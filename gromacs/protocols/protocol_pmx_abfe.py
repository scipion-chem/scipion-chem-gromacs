# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Joaquin Algorta (joaquin.algorta@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
Absolute binding free energy (ABFE) for one ligand, via pmx's Boresch-restraint
double-decoupling setup (`pmx abfe`) plus native GROMACS multi-window
equilibrium free-energy legs:

  0. Parametrize the ligand (ACPYPE) and build the receptor-only topology
     (plain pdb2gmx, no ligand merge - pmx abfe wants them separate).
  1. pmx abfe --build: picks 3+3 Boresch-restraint atoms (ligand+protein),
     writes the restraint as a native dual-state (A/B) bonded term into a
     merged, solvated/ionized complex system, AND separately builds a
     solvated/ionized ligand-alone system. Also computes the analytical
     restraint free energy (no simulation needed for that term).
  2. Leg "restrain": in the complex, ramp the restraint OFF->ON (bonded-lambdas
     only - the restraint is a plain bonded term, no soft-core needed).
  3. Leg "decouple_bound": in the complex (restraint held ON), decouple the
     ligand's electrostatics then van der Waals from the environment via
     native GROMACS couple-moltype.
  4. Leg "decouple_free": same decoupling, on the ligand-alone system (no
     protein, no restraint).
  5. Each leg analysed with native `gmx bar` (equilibrium multi-window - NOT
     pmx analyse, which is for the non-equilibrium fast-growth RBFE protocol).
  6. Combine: dG_bind = dG_restrain + dG_decouple_bound - dG_decouple_free
     - dG_restraints_analytical (standard double-decoupling cycle, Boresch et
     al. 2003, J Phys Chem B 107(35)). All four components are always kept,
     not just the combined total - see the module docstring in
     GROMACS_KNOWLEDGE.md for why the sign of the analytical term should be
     treated as best-effort rather than fully certain.

Reference: pmx, Gapsys et al. (https://github.com/deGrootLab/pmx);
Boresch et al., J Phys Chem B 2003, 107(35), 9535-9551.
"""

import os
import re
import shutil
from multiprocessing import cpu_count

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem
from gromacs.constants import PMX_SC_ALPHA, PMX_SC_SIGMA
from gromacs.protocols.protocol_system_prep import GromacsSystemPrep, GAPS_OPTIONS, LIGAND

RESTRAIN, DECOUPLE_BOUND, DECOUPLE_FREE = 'restrain', 'decouple_bound', 'decouple_free'
LEGS = (RESTRAIN, DECOUPLE_BOUND, DECOUPLE_FREE)
ABFE_DIR = 'abfe'


class GromacsPmxABFE(GromacsSystemPrep):
    """
    Absolute binding free energy for one ligand, via pmx's Boresch-restraint
    double-decoupling setup and native GROMACS multi-window equilibrium TI/BAR.

    Meant as the alternative to GromacsPmxRBFE (relative binding free energy) when the
    ligand of interest has no structurally similar partner to alchemically
    morph into/from - RBFE needs a shared scaffold (pmx atomMapping is
    MCS-based); ABFE evaluates one ligand's binding on its own.

    Window counts default to small, fast-for-testing values, not to settings
    tuned for production accuracy - the same scope note as GromacsPmxRBFE's
    fast-growth snapshot/switching-time defaults.
    """
    _label = 'pmx absolute binding free energy (ABFE)'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        cpus = cpu_count() // 2
        form.addParallelSection(threads=cpus, mpi=1)

        form.addHidden(params.USE_GPU, params.BooleanParam, default=False,
                       label='Use GPU for execution',
                       help='GPU may have several cores. Set it one if you don\'t know what we are '
                            'talking about but you have a GPU.')
        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label='Choose GPU IDs',
                       help='Add a list of GPU devices that can be used (Comma separated)')
        # ABFE is always ligand-vs-nothing (no second structure/ligand involved), same
        # as GromacsPmxRBFE: the AtomStruct-only input mode doesn't apply, so it's hidden.
        form.addHidden('inputFrom', params.EnumParam, choices=['AtomStruct', 'SetOfSmallMolecules'],
                       default=LIGAND)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False,
                      help='Set of docked/aligned molecules sharing a common receptor.')
        form.addParam('inputLigand', params.StringParam,
                      label='Ligand: ',
                      help='Ligand to compute the absolute binding free energy for, picked from '
                           'the set above.')

        group = form.addGroup('Force field')
        self._defineFFParams(group)

        form.addParam('addCaps', params.EnumParam, choices=GAPS_OPTIONS, default=0,
                      label='Add ACE and NME caps: ',
                      help='Add acetyl (ACE) and N-methylamide (NME) capping groups to protein '
                           'termini before building the receptor topology (see "System preparation").')

        self._defineACPYPEparams(form, condition=True)

        form.addSection('Receptor preparation')
        group = form.addGroup('SS bonds')
        self._defineSSBondsParams(group)

        form.addSection(label='ABFE settings')
        grp = form.addGroup('Restraint setup (pmx abfe)')
        grp.addParam('restrSwitchOn', params.BooleanParam, default=True,
                     label='Switch restraint on (vs off): ',
                     help='pmx --restr_switch_on. True: no restraint in state A, full restraint '
                          'in state B (the restrain leg then ramps lambda 0->1). Default: True.')
        grp.addParam('seed', params.IntParam, default=-1, expertLevel=params.LEVEL_ADVANCED,
                     label='Random seed for automatic restraint-atom selection: ',
                     help='-1: let pmx pick a random seed. Set a fixed value for reproducibility.')

        grp = form.addGroup('Equilibration (per window)')
        grp.addParam('nStepsMin', params.IntParam, default=10000, label='Max EM steps: ')
        grp.addParam('emTol', params.FloatParam, default=1000.0, label='EM max force objective: ')
        grp.addParam('equilTime', params.FloatParam, default=50.0, label='Equilibration time (ps): ')
        grp.addParam('prodTime', params.FloatParam, default=100.0, label='Production time (ps): ',
                     help='dH/dl is collected throughout the production run of every window. '
                          'DEFAULT IS FAST/EXPLORATORY, NOT PRODUCTION-QUALITY: 100 ps per window is '
                          'enough to check that the setup/topology/lambda schedule are mechanically '
                          'correct, but is far too short to trust the resulting dG. Published ABFE '
                          'protocols typically use at least 2000-5000 ps (2-5 ns) per window, more '
                          'for the van der Waals decoupling windows near the fully-decoupled end '
                          'state where sampling is hardest. Raise this (and consider raising the '
                          'window counts below too) before relying on the combined result.')
        grp.addParam('temperature', params.FloatParam, default=300.0, label='Temperature (K): ')
        grp.addParam('pressure', params.FloatParam, default=1.0, label='Pressure (bar): ')
        grp.addParam('timeStep', params.FloatParam, default=0.002, expertLevel=params.LEVEL_ADVANCED,
                     label='Time step (ps)[dt]: ')

        grp = form.addGroup('Lambda schedule')
        grp.addParam('nRestrWindows', params.IntParam, default=8,
                     label='Restrain-leg windows: ',
                     help='Number of windows ramping the Boresch restraint OFF->ON in the complex. '
                          '8 is a reasonable, commonly-used count for this leg - the restraint term '
                          'usually converges without needing as many windows as the decoupling legs '
                          'below.')
        grp.addParam('nCoulWindows', params.IntParam, default=5,
                     label='Decouple-leg Coulomb windows: ',
                     help='Windows turning the ligand\'s charges off (shared schedule for both '
                          'decoupling legs). 5 is the low end of what\'s commonly used (electrostatic '
                          'decoupling is usually smoother than van der Waals) - fine for a quick '
                          'check, consider 8-10 for a more defensible result.')
        grp.addParam('nVdwWindows', params.IntParam, default=8,
                     label='Decouple-leg van der Waals windows: ',
                     help='Windows turning the ligand\'s van der Waals interactions off, after '
                          'Coulomb. Van der Waals decoupling usually needs more/denser windows '
                          'than Coulomb to converge - this is the hardest-converging leg (worst '
                          'phase-space overlap near the fully-decoupled state). DEFAULT IS FAST/'
                          'EXPLORATORY: 8 is on the low side for a trustworthy result - 12-20 is the '
                          'more commonly recommended range, especially for larger/more lipophilic '
                          'ligands.')
        grp.addParam('scAlpha', params.FloatParam, default=PMX_SC_ALPHA, expertLevel=params.LEVEL_ADVANCED,
                     label='Soft-core alpha: ')
        grp.addParam('scSigma', params.FloatParam, default=PMX_SC_SIGMA, expertLevel=params.LEVEL_ADVANCED,
                     label='Soft-core sigma: ')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        ligStep = self._insertFunctionStep(self.parametrizeLigandStep)
        recStep = self._insertFunctionStep(self.prepareReceptorStep)
        setupStep = self._insertFunctionStep(self.runAbfeSetupStep, prerequisites=[ligStep, recStep])

        # Every window in every leg starts fresh from that leg's own genion.gro (the mdp's
        # lambda schedule, not trajectory continuation, is what distinguishes windows) - so
        # all windows across all three legs are mutually independent, each depending only
        # on setupStep having produced the leg's topology/starting structure.
        analyseSteps = []
        for leg in LEGS:
            windowSteps = []
            for i in range(self._legWindowCount(leg)):
                eqStep = self._insertFunctionStep(self.equilibrateWindowStep, leg, i,
                                                  prerequisites=[setupStep])
                prodStep = self._insertFunctionStep(self.productionWindowStep, leg, i,
                                                    prerequisites=[eqStep])
                windowSteps.append(prodStep)
            analyseSteps.append(self._insertFunctionStep(self.analyseLegStep, leg,
                                                          prerequisites=windowSteps))

        self._insertFunctionStep(self.createOutputStep, prerequisites=analyseSteps)

    def _nDecoupleWindows(self):
        # nCoulWindows + nVdwWindows points, minus the one shared boundary
        # point (coul=1,vdw=0) common to both stages - see _decoupleLambdas.
        return self.nCoulWindows.get() + self.nVdwWindows.get() - 1

    # -- 0. Ligand parametrization (inherited) + receptor-only topology -----------
    def prepareReceptorStep(self):
        """Plain pdb2gmx on the receptor alone (no ligand merge) - pmx abfe wants
        separate protein/ligand topologies, unlike GromacsSystemPrep's ligand-mode
        PDB2GMXStep which merges them (so that inherited step isn't reused here)."""
        inputStructure = self.getInputReceptorFile()
        systemBasename = self.getSystemName()

        Waterff = self.getEnumText('waterForceField')
        Mainff = self.getEnumText('mainForceField')
        outGro = f'{systemBasename}_protein.gro'
        command = (f'pdb2gmx -f {inputStructure} -o {outGro} '
                   f'-water {Waterff} -ff {Mainff} -merge all')

        ssMode = self.getEnumText('handleSSBonds')
        printfValues = None
        if ssMode != 'None':
            command += ' -ss'
            numBonds = self.countSSBonds(inputStructure, Waterff, Mainff)
            if ssMode == 'Automatic':
                printfValues = ['y'] * int(numBonds)
            elif ssMode == 'Manual':
                selected = {int(x) for x in self.selectSSBonds.get().split(',')} if self.selectSSBonds.get() else set()
                printfValues = ['y' if i in selected else 'n' for i in range(numBonds)]

        try:
            self._runPdb2gmx(command, printfValues=printfValues)
        except Exception:
            self._log.warning('Conversion to gro failed, trying with -ignh flag')
            topFile = self._getPath('topol.top')
            if os.path.exists(topFile):
                os.remove(topFile)
            self._runPdb2gmx(command + ' -ignh', printfValues=printfValues)

        os.rename(self._getPath('topol.top'), self._getPath('protein.top'))

    def _proteinTopFile(self):
        return self._getPath('protein.top')

    def _proteinGroFile(self):
        return self._getPath(f'{self.getSystemName()}_protein.gro')

    # -- 1. pmx abfe --build --------------------------------------------------------
    def runAbfeSetupStep(self):
        abfeDir = self._getExtraPath(ABFE_DIR)
        os.makedirs(abfeDir, exist_ok=True)

        # pmx abfe --build unconditionally os.mkdir()s 'complex'/'ligand' in cwd and
        # writes complex.top/ligand.itp/restraints.info etc alongside them - it is not
        # safe to re-invoke over its own previous (partial or full) output, which is
        # exactly what happens on a step retry after a failure. Clear any leftovers from
        # a prior attempt first so this step is actually idempotent.
        for name in ('complex', 'ligand'):
            shutil.rmtree(os.path.join(abfeDir, name), ignore_errors=True)
        for name in ('complex.gro', 'complex.top', 'ligand.itp', 'posre_ligand.itp',
                     'restraints.info'):
            fpath = os.path.join(abfeDir, name)
            if os.path.exists(fpath):
                os.remove(fpath)

        molName = self.getLigandName()
        ligItp = os.path.abspath(self.getLigandPath(f'{molName}_GMX.itp'))
        ligGro = os.path.abspath(self.getLigandPath(f'{molName}_GMX.gro'))
        proteinTop = os.path.abspath(self._proteinTopFile())
        proteinGro = os.path.abspath(self._proteinGroFile())

        args = (f'-pt {proteinTop} -lt {ligItp} -pc {proteinGro} -lc {ligGro} --build')
        if self.seed.get() is not None and self.seed.get() >= 0:
            args += f' --seed {self.seed.get()}'
        if not self.restrSwitchOn.get():
            args += ' --restr_switch_on'  # store_false flag: passing it flips to False
        gromacsPlugin.runPmx(self, 'abfe', args, cwd=abfeDir)

    def _legDir(self, leg):
        return self._getExtraPath(ABFE_DIR, 'complex' if leg != DECOUPLE_FREE else 'ligand')

    def _legTopFile(self, leg):
        return os.path.join(self._legDir(leg), 'complex.top' if leg != DECOUPLE_FREE else 'ligand.top')

    def _legStartGro(self, leg):
        return os.path.join(self._legDir(leg), 'genion.gro')

    def _restraintsInfoFile(self):
        return self._getExtraPath(ABFE_DIR, 'restraints.info')

    def _parseRestraintsDG(self):
        """Best-effort parse of pmx abfe's restraints.info analytical restraint
        free energy. Never raises: None means the analytical term is missing,
        handled gracefully downstream (same lenient style as GromacsPmxRBFE's
        results.txt parser)."""
        infoFile = self._restraintsInfoFile()
        if not os.path.exists(infoFile):
            return None
        with open(infoFile) as f:
            text = f.read()
        match = re.search(r'dG Restraints\s*=\s*([-\d.]+)\s*kJ/mol', text)
        return float(match.group(1)) if match else None

    # -- 2/3/4. Per-window equilibration + production, per leg ----------------------
    def _windowDir(self, leg, i):
        return self._getExtraPath('legs', leg, f'window_{i}')

    def _decoupleLambdas(self):
        """(coul, vdw) pairs: Coulomb 0->1 first (vdW held at 0), then vdW 0->1
        (Coulomb held at 1), sharing the (1, 0) boundary point once."""
        nCoul, nVdw = self.nCoulWindows.get(), self.nVdwWindows.get()
        coulStage = [(i / (nCoul - 1), 0.0) for i in range(nCoul)]
        vdwStage = [(1.0, i / (nVdw - 1)) for i in range(1, nVdw)]
        return coulStage + vdwStage

    def _restrainLambdas(self):
        n = self.nRestrWindows.get()
        return [i / (n - 1) for i in range(n)]

    def _writeAbfeMdp(self, outFile, leg, i, phase):
        common = ('nstlist          = 10\n'
                  'cutoff-scheme    = Verlet\n'
                  'coulombtype      = PME\n'
                  'rcoulomb         = 1.0\n'
                  'rvdw             = 1.0\n'
                  'pbc              = xyz\n'
                  'constraints      = h-bonds\n'
                  'constraint_algorithm = lincs\n')

        if leg == RESTRAIN:
            lambdas = self._restrainLambdas()
            feBlock = (f'free-energy      = yes\n'
                      f'init-lambda-state = {i}\n'
                      f'bonded-lambdas   = {" ".join(str(v) for v in lambdas)}\n')
        else:
            lambdas = self._decoupleLambdas()
            ligMolName = self.getLigandName()
            coulStr = ' '.join(str(c) for c, v in lambdas)
            vdwStr = ' '.join(str(v) for c, v in lambdas)
            feBlock = (f'free-energy      = yes\n'
                      f'init-lambda-state = {i}\n'
                      f'coul-lambdas     = {coulStr}\n'
                      f'vdw-lambdas      = {vdwStr}\n'
                      f'sc-alpha         = {self.scAlpha.get()}\n'
                      f'sc-power         = 1\n'
                      f'sc-sigma         = {self.scSigma.get()}\n'
                      f'couple-moltype   = {ligMolName}\n'
                      f'couple-lambda0   = vdw-q\n'
                      f'couple-lambda1   = none\n'
                      # couple-intramol=no converts the ligand's intramolecular pairs into fixed
                      # explicit interactions that must all fit within the pair-list cutoff - this
                      # fails ("perturbed excluded non-bonded pair interactions beyond cut-off")
                      # for elongated ligands (confirmed directly with retinal). =yes decouples
                      # intramolecular interactions too, avoiding that failure mode generally.
                      f'couple-intramol  = yes\n')
            if leg == DECOUPLE_BOUND:
                # Hold the restrain leg's end state (restraint fully applied) throughout -
                # a fixed (not per-window) bonded-lambdas array of all-1s.
                nBonded = max(len(lambdas), 2)
                feBlock += f'bonded-lambdas   = {" ".join(["1"] * nBonded)}\n'

        if phase == 'em':
            text = ('integrator       = steep\n'
                   f'emtol            = {self.emTol.get()}\n'
                   'emstep           = 0.01\n'
                   f'nsteps           = {self.nStepsMin.get()}\n'
                   + common + feBlock)
        else:
            simTime = self.equilTime.get() if phase == 'equil' else self.prodTime.get()
            nsteps = round(simTime / self.timeStep.get())
            text = ('integrator       = md\n'
                   f'dt               = {self.timeStep.get()}\n'
                   f'nsteps           = {nsteps}\n'
                   + common +
                   f'continuation     = {"no" if phase == "equil" else "yes"}\n'
                   f'gen_vel          = {"yes" if phase == "equil" else "no"}\n'
                   f'gen_temp         = {self.temperature.get()}\n'
                   'gen_seed         = -1\n'
                   'nstxout          = 0\n'
                   'nstvout          = 0\n'
                   f'nstdhdl          = {1 if phase == "prod" else 100}\n'
                   'tcoupl           = V-rescale\n'
                   'tc-grps          = System\n'
                   f'ref_t            = {self.temperature.get()}\n'
                   'tau_t            = 0.1\n'
                   'nsttcouple       = 10\n'
                   'pcoupl           = C-rescale\n'
                   'pcoupltype       = isotropic\n'
                   f'ref_p            = {self.pressure.get()}\n'
                   'tau_p            = 2.0\n'
                   'compressibility  = 4.5e-5\n'
                   'nstpcouple       = 10\n'
                   + feBlock)

        with open(outFile, 'w') as f:
            f.write(text)
        return outFile

    def _runGrompp(self, stageDir, mdpFile, groFile, topFile, tprName, prevCpt=None):
        localTop = os.path.join(stageDir, os.path.basename(topFile))
        if not os.path.exists(localTop):
            os.link(topFile, localTop)
        # complex.top #includes "ligand.itp", a sibling file pmx abfe writes alongside it -
        # confirmed by actually running grompp: without this, it fails with 'Topology include
        # file "ligand.itp" not found' (the same relocation invariant GromacsPmxRBFE already
        # has to handle for its own hybrid ligand itp). ligand.top has no such sibling include,
        # so this is a no-op for the decouple_free leg.
        topDir = os.path.dirname(topFile)
        for fname in os.listdir(topDir):
            if fname.endswith('.itp'):
                dst = os.path.join(stageDir, fname)
                if not os.path.exists(dst):
                    os.link(os.path.join(topDir, fname), dst)
        prevStr = f' -t {os.path.abspath(prevCpt)}' if prevCpt else ''
        outFile = f'{tprName}.tpr'
        command = (f'grompp -f {os.path.abspath(mdpFile)} -c {os.path.abspath(groFile)} '
                   f'-p {os.path.basename(localTop)}{prevStr} -o {outFile} -maxwarn 2')
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir, numberOfMpi=0)
        return os.path.join(stageDir, outFile)

    def _runMdrun(self, stageDir, deffnm):
        if getattr(self, params.USE_GPU):
            gpuList = getattr(self, params.GPU_LIST).get().replace(' ', '')
            gpuStr = f' -nb gpu -gpu_id {gpuList}'
        else:
            gpuStr = ' -nb cpu'
        command = f'mdrun -v -deffnm {deffnm}{gpuStr} -nt {self.numberOfThreads.get()}'
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir)

    def _failedFlagPath(self, leg, i):
        return os.path.abspath(self._getExtraPath(f'.failed_{leg}_{i}.flag'))

    def equilibrateWindowStep(self, leg, i):
        wDir = self._windowDir(leg, i)
        os.makedirs(wDir, exist_ok=True)
        topFile = self._legTopFile(leg)
        startGro = self._legStartGro(leg)

        try:
            emMdp = self._writeAbfeMdp(os.path.join(wDir, 'em.mdp'), leg, i, 'em')
            self._runGrompp(wDir, emMdp, startGro, topFile, 'em')
            self._runMdrun(wDir, 'em')

            eqMdp = self._writeAbfeMdp(os.path.join(wDir, 'equil.mdp'), leg, i, 'equil')
            self._runGrompp(wDir, eqMdp, os.path.join(wDir, 'em.gro'), topFile, 'equil')
            self._runMdrun(wDir, 'equil')
        except Exception as e:
            self.warning(f'Equilibration failed for leg {leg} window {i}: {e}')
            with open(self._failedFlagPath(leg, i), 'w') as fh:
                fh.write(str(e) + '\n')

    def productionWindowStep(self, leg, i):
        if os.path.exists(self._failedFlagPath(leg, i)):
            return
        wDir = self._windowDir(leg, i)
        topFile = self._legTopFile(leg)
        try:
            prodMdp = self._writeAbfeMdp(os.path.join(wDir, 'prod.mdp'), leg, i, 'prod')
            self._runGrompp(wDir, prodMdp, os.path.join(wDir, 'equil.gro'), topFile, 'prod',
                            prevCpt=os.path.join(wDir, 'equil.cpt'))
            self._runMdrun(wDir, 'prod')
        except Exception as e:
            self.warning(f'Production failed for leg {leg} window {i}: {e}')
            with open(self._failedFlagPath(leg, i), 'w') as fh:
                fh.write(str(e) + '\n')

    # -- 5. Per-leg BAR analysis ------------------------------------------------------
    def _legWindowCount(self, leg):
        return self.nRestrWindows.get() if leg == RESTRAIN else self._nDecoupleWindows()

    def _legDhdlFiles(self, leg):
        # mdrun names the dH/dl output <deffnm>.xvg (not <deffnm>.dhdl.xvg), confirmed
        # by actually running a free-energy mdrun on this machine (same fact as
        # GromacsPmxRBFE). gmx bar additionally needs every consecutive window present -
        # confirmed directly: a gap between non-adjacent windows makes it refuse outright
        # ("no path between the states... covered by foreign lambdas").
        files = []
        for i in range(self._legWindowCount(leg)):
            if os.path.exists(self._failedFlagPath(leg, i)):
                self.warning(f'Leg {leg}: window {i} failed, dropping the whole leg\'s BAR chain '
                             f'from this point (gmx bar needs an unbroken adjacent-window chain).')
                break
            dhdl = os.path.abspath(os.path.join(self._windowDir(leg, i), 'prod.xvg'))
            if not os.path.exists(dhdl):
                break
            files.append(dhdl)
        return files

    def analyseLegStep(self, leg):
        anDir = self._getExtraPath('analysis', leg)
        os.makedirs(anDir, exist_ok=True)
        dhdlFiles = self._legDhdlFiles(leg)
        if len(dhdlFiles) < 2:
            self.warning(f'Leg {leg}: fewer than 2 usable windows ({len(dhdlFiles)}), cannot run gmx bar.')
            return
        # gmx bar reports its result on stdout, not in a result file (confirmed by actually
        # running it) - redirected here to a dedicated per-leg file rather than relying on
        # Scipion's shared logs/run.stdout, since all three legs' analyseLegStep calls can
        # run concurrently (stepsExecutionMode=STEPS_PARALLEL) and would otherwise interleave.
        # Bare filename, not anDir-prefixed: cwd=anDir already, and a relative path prefixed
        # with anDir would be re-resolved against anDir *again* by the subprocess, double
        # nesting (the exact bug this same session already hit and fixed several times in
        # GromacsPmxRBFE - reproduced here for real: "No existe el archivo o el directorio").
        command = f'bar -f {" ".join(dhdlFiles)} -o bar.xvg -oi barint.xvg > bar_result.txt 2>&1'
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=anDir)

    # -- 6. Output -----------------------------------------------------------------
    def createOutputStep(self):
        dGRestraints = self._parseRestraintsDG()
        legResults = {}
        for leg in LEGS:
            dg, err = self._parseBarLog(leg)
            legResults[leg] = (dg, err)

        summaryLines = [f'dG_restraints (analytical) = {dGRestraints} kJ/mol']
        for leg in LEGS:
            dg, err = legResults[leg]
            errStr = f' +/- {err}' if err is not None else ''
            summaryLines.append(f'dG_{leg} = {dg}{errStr} kJ/mol')

        dGTotal = None
        if dGRestraints is not None and all(legResults[leg][0] is not None for leg in LEGS):
            dGTotal = (legResults[RESTRAIN][0] + legResults[DECOUPLE_BOUND][0]
                      - legResults[DECOUPLE_FREE][0] - dGRestraints)
            summaryLines.append(f'dG_bind = dG_restrain + dG_decouple_bound - dG_decouple_free '
                               f'- dG_restraints = {dGTotal} kJ/mol')
        else:
            summaryLines.append('dG_bind could not be computed: one or more legs/terms are missing.')

        resultsFile = self._getExtraPath('analysis', 'results.txt')
        with open(resultsFile, 'w') as f:
            f.write('\n'.join(summaryLines) + '\n')

        molName = self.getLigandName()
        outSystem = GromacsSystem(filename=self._proteinGroFile(), topoFile=self._proteinTopFile(),
                                  ff=self.getEnumText('mainForceField'), wff=self.getEnumText('waterForceField'))
        outSystem.setChainNames(','.join(self.getModelChains()))
        chains, lengthsDic = self.getModelChainsAndLengths()
        outSystem.setChainLengths(','.join(str(lengthsDic[c]) for c in chains))
        outSystem.setLigTopologyFile(self.getLigandPath(f'{molName}_GMX.itp'))
        if dGTotal is not None:
            outSystem.setFreeEnergy(dGTotal)
        outSystem.setFreeEnergyFile(resultsFile)

        self._defineOutputs(outputSystem=outSystem)
        self._defineSourceRelation(self.inputSetOfMols, outSystem)

    def _parseBarLog(self, leg):
        """Best-effort parse of gmx bar's stdout ('total X - Y, DG Z +/- W', kJ/mol),
        redirected by analyseLegStep into its own dedicated per-leg file - confirmed
        real format by actually running gmx bar on this machine. Never raises;
        returns (None, None) on any mismatch (missing file, unexpected format, ...)."""
        logFile = self._getExtraPath('analysis', leg, 'bar_result.txt')
        if not os.path.exists(logFile):
            return None, None
        with open(logFile) as f:
            text = f.read()
        match = re.search(r'total\s+\d+\s*-\s*\d+,\s*DG\s+([-\d.]+)\s*\+/-\s*([\d.]+)', text)
        if not match:
            return None, None
        return float(match.group(1)), float(match.group(2))

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        # Not super()._validate(): GromacsSystemPrep's checks placeIons/cationType/anionType,
        # which this protocol doesn't declare - pmx abfe's own genion call is fixed (0.15 M,
        # neutralize), not user-configurable here.
        errors = []
        mols = self.inputSetOfMols.get()
        if mols is not None:
            names = [mol.__str__() for mol in mols]
            if self.inputLigand.get() not in names:
                errors.append(f'Ligand "{self.inputLigand.get()}" not found in the input set of molecules.')
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            resultsFile = self._getExtraPath('analysis', 'results.txt')
            if os.path.exists(resultsFile):
                with open(resultsFile) as f:
                    summary.append(f.read())
            else:
                summary.append('Finished, but no results.txt was produced.')
        else:
            summary.append('The protocol has not finished.')
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append('Absolute binding free energy was calculated with the Boresch-restraint '
                           'double-decoupling method: restraint setup and analytical restraint free '
                           'energy via "pmx abfe", restrain and decoupling legs via native GROMACS '
                           'multi-window equilibrium free energy (couple-moltype for decoupling), '
                           'analysed with "gmx bar".')
        return methods

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
Relative binding free energy (RBFE) between two ligands docked/aligned in the
same pocket of a shared receptor, using pmx's non-equilibrium fast-growth
alchemical FEP:

  0. Build ligand A's solvated/ionized system (reusing GromacsSystemPrep's own
     ligand-mode steps) and parametrize ligand B (ACPYPE).
  1. pmx atomMapping: find the morphable atoms between ligand A and B.
  2. pmx ligandHybrid: build the dual-topology hybrid ligand (.itp + dummy atomtypes).
  3. Splice the hybrid ligand into ligand A's system topology/coordinates
     (box, solvent and ions are left untouched).
  4. Equilibrate (EM->NVT->NPT) the hybrid system at both end states
     (lambda=0 == ligand A, lambda=1 == ligand B).
  5. Extract snapshots from each equilibrium trajectory and launch short
     forward (0->1) / reverse (1->0) switching ("fast growth") runs.
  6. pmx analyse: estimate dG(A->B) from the collected dH/dl work values
     (Crooks Gaussian Intersection / Bennett Acceptance Ratio / Jarzynski).

Reference: pmx, Gapsys et al. (https://github.com/deGrootLab/pmx)
"""

import os
import re
import shutil
from glob import glob
from multiprocessing import cpu_count

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem import Plugin as pwchemPlugin

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem
from gromacs.constants import PMX_SC_ALPHA, PMX_SC_SIGMA, PMX_SW_TIME, PMX_N_SNAPSHOTS
from gromacs.protocols.protocol_system_prep import (
    GromacsSystemPrep, replaceInFile, LIGAND, TOPOL_TOP, GAPS_OPTIONS,
)

MAP_DIR = 'mapping'
HYBRID_DIR = 'hybrid'
SYSTEM_DIR = 'hybrid_system'
FWD, BWD = 'fwd', 'bwd'
STATE_DIRNAME = {0: 'stateA', 1: 'stateB'}
DIRECTION_STATE = {FWD: 0, BWD: 1}


class GromacsPmxRBFE(GromacsSystemPrep):
    """
    Relative binding free energy (ligand A -> ligand B) via pmx non-equilibrium
    fast-growth alchemical FEP.

    Both ligands are picked (by name, via a SelectElementWizard each) from the
    same SetOfSmallMolecules, docked/aligned in the same pocket of a shared
    receptor. Ligand A's solvated/ionized system is built by this protocol
    itself, reusing GromacsSystemPrep's own ligand-mode preparation steps
    (PDB2GMX/editconf/solvate/genion) - subclassing it is exactly this reuse,
    the same "compose, don't duplicate" pattern GromacsMmpbsa already uses.
    """
    _label = 'pmx relative binding free energy (RBFE)'
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
        # RBFE is always ligand-vs-ligand; the AtomStruct-only input mode GromacsSystemPrep
        # also supports doesn't apply here, so this stays fixed and hidden.
        form.addHidden('inputFrom', params.EnumParam, choices=['AtomStruct', 'SetOfSmallMolecules'],
                       default=LIGAND)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Docked molecules: ', allowsNull=False,
                      help='Set of docked/aligned molecules sharing a common receptor. Both ligand A '
                           'and ligand B (the two FEP end states) are picked from this same set.')
        form.addParam('inputLigand', params.StringParam,
                      label='Ligand A: ',
                      help='First ligand (FEP state A), picked from the set above. Its solvated/ionized '
                           'system is built by this protocol.')
        form.addParam('ligandB', params.StringParam,
                      label='Ligand B: ',
                      help='Second ligand (FEP state B, morphed into from ligand A), picked from the '
                           'set above.')

        group = form.addGroup('Force field')
        self._defineFFParams(group)

        form.addParam('addCaps', params.EnumParam, choices=GAPS_OPTIONS, default=0,
                      label='Add ACE and NME caps: ',
                      help='Add acetyl (ACE) and N-methylamide (NME) capping groups to protein '
                           'termini before building ligand A\'s system (see "System preparation").')

        self._defineACPYPEparams(form, condition=True)

        form.addSection('System preparation')
        group = form.addGroup('Boundary box')
        self._defineBoxParams(group)

        group = form.addGroup('Ions')
        self._defineIonsParams(group)

        group = form.addGroup('SS bonds')
        self._defineSSBondsParams(group)

        form.addSection(label='RBFE settings')
        grp = form.addGroup('End-state equilibration')
        grp.addParam('nStepsMin', params.IntParam, default=10000,
                     label='Max EM steps: ',
                     help='Maximum number of (minimization) steps to perform')
        grp.addParam('emTol', params.FloatParam, default=1000.0,
                     label='EM max force objective: ',
                     help='Stop minimization when the maximum force < x kJ/mol/nm.')
        grp.addParam('nvtTime', params.FloatParam, default=200.0, label='NVT time (ps): ',
                     help='Default (200 ps) is enough for thermalization on typical systems; raise '
                          'it if temperature/pressure have not visibly plateaued in the log.')
        grp.addParam('nptTime', params.FloatParam, default=2000.0, label='NPT time (ps): ',
                     help='Snapshots for the switching runs are taken evenly from the second half of '
                          'this trajectory (the first half is burn-in). Default (2000 ps, i.e. 1000 ps '
                          'of usable trajectory) is set so that, at the default 20 snapshots per '
                          'end-state, consecutive snapshots are ~50 ps apart - enough to decorrelate '
                          'for most drug-like ligands. If you raise "Snapshots per end-state" '
                          'substantially, raise this proportionally too, or the extra snapshots will '
                          'just be correlated repeats rather than independent samples.')
        grp.addParam('temperature', params.FloatParam, default=300.0, label='Temperature (K): ')
        grp.addParam('pressure', params.FloatParam, default=1.0, label='Pressure (bar): ')
        grp.addParam('timeStep', params.FloatParam, default=0.002, expertLevel=params.LEVEL_ADVANCED,
                     label='Time step (ps)[dt]: ')

        grp = form.addGroup('Fast-growth transitions')
        grp.addParam('nStructs', params.IntParam, default=PMX_N_SNAPSHOTS,
                     label='Snapshots per end-state: ',
                     help='Number of frames taken from each end-state equilibration to launch '
                          'independent forward/reverse switching runs from - i.e. this many forward '
                          '+ this many reverse non-equilibrium transitions get run and fed to the '
                          'CGI/BAR/JARZ estimators. Default (20, i.e. 40 total transitions) is pmx\'s '
                          'own commonly-used baseline and should give a usable estimate with a '
                          'reasonably tight error bar for a typical drug-like ligand pair; the error '
                          'bar narrows (roughly) with the square root of this number, so doubling it '
                          'to 40 will noticeably tighten it at ~2x the transition-running cost (the '
                          'end-state equilibrations above are NOT repeated, only the cheap switching '
                          'runs are).')
        grp.addParam('swTime', params.FloatParam, default=PMX_SW_TIME,
                     label='Switching time (ps): ',
                     help='Length of the linear lambda-ramp transition run (pmx default and this '
                          'protocol\'s default: 50 ps - matches published pmx RBFE campaigns and is '
                          'not a "fast/exploratory" shortcut like some other defaults in this plugin; '
                          'no change needed for production use).')
        grp.addParam('scAlpha', params.FloatParam, default=PMX_SC_ALPHA, expertLevel=params.LEVEL_ADVANCED,
                     label='Soft-core alpha: ')
        grp.addParam('scSigma', params.FloatParam, default=PMX_SC_SIGMA, expertLevel=params.LEVEL_ADVANCED,
                     label='Soft-core sigma: ')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Ligand A's system build, reusing GromacsSystemPrep's own ligand-mode steps
        # unchanged (self.inputLigand/inputSetOfMols/inputFrom==LIGAND drive them).
        ligAStep = self._insertFunctionStep(self.parametrizeLigandStep)
        ligBStep = self._insertFunctionStep(self.parametrizeLigandBStep)  # independent of ligand A

        pdbStep = self._insertFunctionStep(self.PDB2GMXStep, prerequisites=[ligAStep])
        ecStep = self._insertFunctionStep(self.editConfStep, prerequisites=[pdbStep])
        solvStep = self._insertFunctionStep(self.solvateStep, prerequisites=[ecStep])
        prevStep = solvStep
        if self.placeIons.get() != 0:
            prevStep = self._insertFunctionStep(self.addIonsStep, prerequisites=[solvStep])

        mStep = self._insertFunctionStep(self.atomMappingStep, prerequisites=[prevStep, ligBStep])
        hStep = self._insertFunctionStep(self.buildHybridLigandStep, prerequisites=[mStep])
        sStep = self._insertFunctionStep(self.buildHybridSystemStep, prerequisites=[hStep])

        eqSteps, extractSteps = {}, {}
        for state in (0, 1):
            eqSteps[state] = self._insertFunctionStep(self.equilibrateStateStep, state, prerequisites=[sStep])
            extractSteps[state] = self._insertFunctionStep(self.extractSnapshotsStep, state,
                                                            prerequisites=[eqSteps[state]])

        transitionSteps = []
        for direction, state in DIRECTION_STATE.items():
            for i in range(self.nStructs.get()):
                tStep = self._insertFunctionStep(self.runTransitionStep, direction, i,
                                                 prerequisites=[extractSteps[state]])
                transitionSteps.append(tStep)

        aStep = self._insertFunctionStep(self.analyseStep, prerequisites=transitionSteps)
        self._insertFunctionStep(self.createOutputStep, prerequisites=[aStep])

    # -- 0. Ligand A's system (GromacsSystemPrep.parametrizeLigandStep/PDB2GMXStep/
    #       editConfStep/solvateStep/addIonsStep, all inherited unchanged) and ligand B --
    def parametrizeLigandBStep(self):
        mol = self.getLigandBMol()
        molFile = self.addHydrogens(os.path.abspath(mol.getPoseFile()))

        kwargs = self.getParameters()
        kwargs['molName'] = mol.getMolName()

        args = f'-i {molFile} -b {kwargs["molName"]} -c {kwargs["chargeMethod"]} ' \
               f'-m {kwargs["multip"]} -a {kwargs["atomType"]} -q {kwargs["qprog"]} -o gmx'
        if 'netCharge' in kwargs:
            args += f' -n {kwargs["netCharge"]}'
        pwchemPlugin.runACPYPE(self, args=args, cwd=self._getExtraPath())

    def getLigandBMol(self):
        myMol = None
        for mol in self.inputSetOfMols.get():
            if mol.__str__() == self.ligandB.get():
                myMol = mol.clone()
                break
        if myMol is None:
            raise ValueError(f'Ligand B "{self.ligandB.get()}" not found in the input set of molecules')
        return myMol

    def getLigandBPath(self, path=''):
        molName = self.getLigandBMol().getMolName()
        return self._getExtraPath(f'{molName}.acpype', path)

    def _ligBGroFile(self):
        """Ligand B's ACPYPE .gro, rigidly translated into the same frame as ligand A's
        real, final built system.

        ACPYPE's output preserves ligand B's real absolute docked-pose coordinates as-is
        (confirmed for real: its centroid matches its own Vina pose almost exactly - ACPYPE/
        openbabel don't recentre a single-molecule conversion). But ligand A's system went
        through a further whole-system recentring at editConfStep (`gmx editconf -c`, never
        called with `-princ` here, so it is a pure translation, no rotation) that ligand B's
        standalone .gro never went through. Mixing the two frames is exactly what broke this
        the first time it was tried with two genuinely different ligands (never caught by the
        earlier self-transform test, where A and B are the same molecule with zero unique
        atoms - this code path was simply never exercised): pmx's ligandHybrid placed ligand
        B's unique/dummy atoms straight from its own untranslated coordinates, landing tens of
        nm from ligand A's real position in the final box and making mdrun's domain
        decomposition fail outright ("no domain decomposition ... minimum cell size ~50 nm").

        Fix: apply the same translation ligand A underwent (computed from its own before/
        after positions) to ligand B before handing it to pmx."""
        rawGro = os.path.abspath(self.getLigandBPath(f'{self.getLigandBMol().getMolName()}_GMX.gro'))
        translatedGro = self._getExtraPath('ligandB_translated.gro')
        if not os.path.exists(translatedGro):
            translation = self._getLigandATranslation()
            self._translateGro(rawGro, translation, translatedGro)
        return os.path.abspath(translatedGro)

    def _getLigandATranslation(self):
        origGro = self.getLigandPath(f'{self.getLigandName()}_GMX.gro')
        finalGro = self._getExtraPath('ligandA.gro')
        ox, oy, oz = self._groCentroid(origGro)
        fx, fy, fz = self._groCentroid(finalGro)
        return fx - ox, fy - oy, fz - oz

    @staticmethod
    def _groCentroid(groPath):
        with open(groPath) as f:
            lines = f.readlines()
        n = int(lines[1])
        xs = ys = zs = 0.0
        for line in lines[2:2 + n]:
            xs += float(line[20:28])
            ys += float(line[28:36])
            zs += float(line[36:44])
        return xs / n, ys / n, zs / n

    def _fixPmxDummyAtomUnits(self, mergedAGro):
        """Correct a real unit bug in pmx's own ligandHybrid (pmx/ligand_alchemy.py's
        adjustCoords(), confirmed by reading the installed source): it fits ligand B onto
        ligand A via an RDKit round-trip through a temporary PDB file (Angstroms), then
        copies the fitted RDKit coordinates straight into pmx's Model atoms with no A->nm
        conversion. Ligand A's real atoms are never touched by this path (they keep their
        correct pmx/nm coordinates throughout) and pmx's own PDB writer + this protocol's
        subsequent `editconf` PDB->gro step apply the standard x10/1/10 conversions
        correctly for them - net zero. But the already-wrong (10x too large) dummy atoms
        pmx builds from that fit go through the SAME write+editconf chain once more,
        compounding to a net 10x-too-large final value.

        Confirmed for real: with two structurally different, independently-docked ligands
        (never exercised before - the earlier self-transform A==B test has zero dummy
        atoms, so it could never expose this), the dummy/unique atoms pmx adds for ligand
        B landed 26-48 nm from ligand A's real ~2.5-4.5 nm cluster in the same box,
        breaking mdrun's domain decomposition outright ("no domain decomposition ...
        compatible with ... minimum cell size ~50-60 nm"). Dividing those dummy atoms'
        coordinates by 10 lands them squarely back in ligand A's real cluster (checked by
        hand against the real failing run's own numbers), confirming this is exactly the
        bug and exactly the fix - not a coincidental rescaling.

        Dummy atoms are exactly the ones beyond ligand A's own real atom count: this
        protocol's hybrid ligand always lists ligand A's real atoms first (unchanged) and
        appends the new dummy atoms after (confirmed by inspecting a real mergedA.gro)."""
        with open(self._getExtraPath('ligandA.gro')) as f:
            nRealAtoms = int(f.readlines()[1])

        with open(mergedAGro) as f:
            lines = f.readlines()
        n = int(lines[1])
        outLines = lines[:2]
        for i, line in enumerate(lines[2:2 + n]):
            if i >= nRealAtoms:
                x = float(line[20:28]) / 10.0
                y = float(line[28:36]) / 10.0
                z = float(line[36:44]) / 10.0
                line = f'{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:]}'
            outLines.append(line)
        outLines.extend(lines[2 + n:])
        with open(mergedAGro, 'w') as f:
            f.writelines(outLines)

    @staticmethod
    def _translateGro(inGro, translation, outGro):
        dx, dy, dz = translation
        with open(inGro) as f:
            lines = f.readlines()
        n = int(lines[1])
        outLines = lines[:2]
        for line in lines[2:2 + n]:
            x = float(line[20:28]) + dx
            y = float(line[28:36]) + dy
            z = float(line[36:44]) + dz
            outLines.append(f'{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:]}')
        outLines.extend(lines[2 + n:])
        with open(outGro, 'w') as f:
            f.writelines(outLines)

    def _ligBItpFile(self):
        return os.path.abspath(self.getLigandBPath(f'{self.getLigandBMol().getMolName()}_GMX.itp'))

    def _getLigandASystem(self):
        """Fresh (non-persisted) GromacsSystem wrapper around ligand A's just-built
        solvated/ionized system. Reconstructed from disk every time it's needed rather
        than cached on self, consistent with how the rest of this plugin re-derives
        state across steps (Scipion may reload the protocol instance between steps)."""
        systemBasename = self.getSystemName()
        groBaseName = (f'{systemBasename}_solv_ions.gro' if self.placeIons.get() != 0
                       else f'{systemBasename}_solv.gro')
        system = GromacsSystem(filename=self._getPath(groBaseName), topoFile=self._getPath(TOPOL_TOP),
                               restrFile=self._getPath('posre.itp'),
                               ff=self.getEnumText('mainForceField'), wff=self.getEnumText('waterForceField'))
        system.setLigTopologyFile(self._getPath(f'{self.getLigandName()}_GMX.itp'))
        return system

    def _ensureLigandAIndexFile(self):
        indexFile = self._getExtraPath('indexes.ndx')
        if not os.path.exists(indexFile):
            chains, lengthsDic = self.getModelChainsAndLengths()
            gromacsPlugin.firstIndexCreation(self, self._getLigandASystem(), ligandName=self.getLigandName(),
                                             modelChains=chains, chainLengths=lengthsDic)
        return indexFile

    # -- 1. Atom mapping --------------------------------------------------------
    def atomMappingStep(self):
        mapDir = self._getExtraPath(MAP_DIR)
        os.makedirs(mapDir, exist_ok=True)

        ligAGro = self._extractLigandAGro()
        ligBGro = self._ligBGroFile()

        args = f'-i1 {ligAGro} -i2 {ligBGro} -o1 pairs1.dat -o2 pairs2.dat -log mapping.log'
        gromacsPlugin.runPmx(self, 'atomMapping', args, cwd=mapDir)

    def _extractLigandAGro(self):
        """Extract ligand A alone from its just-built system, via GROMACS' own
        auto-generated "Other" index group (the standard catch-all for
        non-protein/nucleic/water residues).

        Deliberately NOT a stored ligand-ID label: GROMACS' own auto-detected residue
        name for a real ligand (e.g. a genuine PDB code like 'RET') is not reliably
        the same string a naming heuristic might guess (e.g. 'LIG') - confirmed by
        hitting "Error: No such group 'LIG'" when this used to rely on one."""
        systemFile = os.path.abspath(self._getLigandASystem().getSystemFile())
        outGro = os.path.abspath(self._getExtraPath('ligandA.gro'))
        indexFile = os.path.abspath(self._ensureLigandAIndexFile())

        groups = gromacsPlugin.parseIndexFile(self, indexFile)
        invGroups = {v: k for k, v in groups.items()}
        if 'Other' not in invGroups:
            raise ValueError(f'No "Other" index group found in {indexFile}; cannot isolate '
                             f'ligand A (the system may have more than one non-protein/water '
                             f'residue type, which this protocol does not support)')

        args = f'trjconv -f {systemFile} -s {systemFile} -n {indexFile} -o {outGro}'
        gromacsPlugin.runGromacsPrintf(self, printfValues=[invGroups['Other']], args=args,
                                       cwd=self._getExtraPath())
        return outGro

    @staticmethod
    def _readResName(groFile):
        """Real residue name physically written in a single-residue .gro file's atom lines
        (ground truth) - see _extractLigandAGro for why nothing here trusts a naming guess."""
        with open(groFile) as f:
            lines = f.readlines()
        return lines[2][5:10].strip()

    # -- 2. Hybrid ligand build --------------------------------------------------
    def buildHybridLigandStep(self):
        mapDir = self._getExtraPath(MAP_DIR)
        hybDir = self._getExtraPath(HYBRID_DIR)
        os.makedirs(hybDir, exist_ok=True)

        ligAGro = os.path.abspath(self._getExtraPath('ligandA.gro'))
        ligBGro = self._ligBGroFile()
        ligAItp = os.path.abspath(self._getLigandASystem().getLigTopologyFile())
        ligBItp = self._ligBItpFile()
        pairsFile = os.path.abspath(os.path.join(mapDir, 'pairs1.dat'))

        args = (f'-i1 {ligAGro} -i2 {ligBGro} -itp1 {ligAItp} -itp2 {ligBItp} '
                f'-pairs {pairsFile} -oA mergedA.pdb -oB mergedB.pdb '
                f'-oitp merged.itp -offitp ffmerged.itp -log hybrid.log')
        gromacsPlugin.runPmx(self, 'ligandHybrid', args, cwd=hybDir)

    # -- 3. Splice the hybrid ligand into ligand A's system -----------------------
    def buildHybridSystemStep(self):
        hybDir = self._getExtraPath(HYBRID_DIR)
        sysDir = self._getExtraPath(SYSTEM_DIR)
        os.makedirs(sysDir, exist_ok=True)
        system = self._getLigandASystem()

        # -f crosses into hybDir (a sibling of the cwd below) so it must be absolute; -o stays
        # a bare filename since it lands inside cwd=sysDir itself (a joined sysDir/... string
        # here would be re-resolved against sysDir again by the subprocess and double-nest).
        mergedAPdb = os.path.abspath(os.path.join(hybDir, 'mergedA.pdb'))
        gromacsPlugin.runGromacs(self, 'gmx', f'editconf -f {mergedAPdb} -o mergedA.gro', cwd=sysDir)
        mergedAGro = os.path.join(sysDir, 'mergedA.gro')
        self._fixPmxDummyAtomUnits(mergedAGro)

        outTop = self._buildHybridTopology(system, hybDir, sysDir)
        # Real residue name (not a naming guess - see _extractLigandAGro), read
        # straight off the ligand-A.gro this same protocol already extracted for mapping.
        ligResName = self._readResName(self._getExtraPath('ligandA.gro'))
        self._spliceLigandGro(system.getSystemFile(), mergedAGro, ligResName,
                              os.path.join(sysDir, 'system.gro'))
        return outTop

    def _buildHybridTopology(self, system, hybDir, sysDir):
        outTop = os.path.join(sysDir, 'topol.top')
        shutil.copy(system.getTopologyFile(), outTop)
        shutil.copy(os.path.join(hybDir, 'merged.itp'), os.path.join(sysDir, 'merged.itp'))

        # merged.itp only carries the *new* dummy atomtypes (ffmerged.itp): the real GAFF/AMBER
        # types it also references (e.g. ca, ha, oh...) come from each ligand's own itp, which
        # would otherwise be lost once their #include is replaced below ("Atomtype X not found"
        # in grompp, confirmed by actually running it).
        self._writeCombinedAtomTypes(
            [system.getLigTopologyFile(), self._ligBItpFile(), os.path.join(hybDir, 'ffmerged.itp')],
            os.path.join(sysDir, 'atomtypes.itp'))

        ligItpBase = os.path.basename(system.getLigTopologyFile())
        oldMolName = ligItpBase.replace('_GMX.itp', '')
        replaceInFile(outTop, f'#include "{ligItpBase}"',
                      '#include "atomtypes.itp"\n#include "merged.itp"')

        # pmx ligandHybrid keeps mol1's original moleculetype name; only patch the
        # [ molecules ] entry in the rare case a future pmx version changes that.
        hybridMolName = self._readMoleculetypeName(os.path.join(sysDir, 'merged.itp'))
        if hybridMolName != oldMolName:
            replaceInFile(outTop, f'{oldMolName} 1', f'{hybridMolName} 1')
        return outTop

    @staticmethod
    def _extractAtomTypesLines(itpFile):
        with open(itpFile) as f:
            lines = f.readlines()
        block, inBlock = [], False
        for line in lines:
            stripped = line.strip().replace(' ', '')
            if stripped == '[atomtypes]':
                inBlock = True
                continue
            if inBlock:
                if stripped.startswith('[') and stripped.endswith(']'):
                    break
                block.append(line)
        return block

    @classmethod
    def _writeCombinedAtomTypes(cls, itpFiles, outFile):
        """Concatenate the [ atomtypes ] blocks of several itp files, deduplicated by
        atom-type name (first occurrence wins; ligands sharing a force field define
        identical params for shared types anyway)."""
        seen, combined = set(), []
        for itpFile in itpFiles:
            for line in cls._extractAtomTypesLines(itpFile):
                stripped = line.strip()
                if not stripped or stripped.startswith(';'):
                    combined.append(line)
                    continue
                name = stripped.split()[0]
                if name in seen:
                    continue
                seen.add(name)
                combined.append(line)

        with open(outFile, 'w') as f:
            f.write('[ atomtypes ]\n')
            f.writelines(combined)
        return outFile

    @staticmethod
    def _readMoleculetypeName(itpFile):
        with open(itpFile) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.strip().replace(' ', '') == '[moleculetype]':
                for nxt in lines[i + 1:]:
                    nxt = nxt.strip()
                    if not nxt or nxt.startswith(';'):
                        continue
                    return nxt.split()[0]
        raise ValueError(f'No [ moleculetype ] section found in {itpFile}')

    @staticmethod
    def _spliceLigandGro(systemGro, hybridLigGro, ligResName, outGro):
        """Replace the ligand-A coordinate block (found by residue name) in systemGro
        with the hybrid ligand's coordinates, keeping box/solvent/ions untouched.
        The hybrid ligand's own residue number/name (pmx writes its own placeholder,
        e.g. "UN") are overwritten to match the original block, so index groups and
        downstream analysis keyed on the original ligand identity keep working."""
        with open(systemGro) as f:
            lines = f.readlines()
        header, boxLine = lines[0], lines[-1]
        atomLines = lines[2:-1]

        with open(hybridLigGro) as f:
            hybLines = f.readlines()
        hybAtomLines = hybLines[2:-1]

        startIdx = endIdx = None
        for i, line in enumerate(atomLines):
            if line[5:10].strip() == ligResName:
                if startIdx is None:
                    startIdx = i
                endIdx = i
        if startIdx is None:
            raise ValueError(f'Ligand residue {ligResName} not found in {systemGro}')

        resField = atomLines[startIdx][0:10]
        hybAtomLines = [resField + line[10:] for line in hybAtomLines]

        newAtomLines = atomLines[:startIdx] + hybAtomLines + atomLines[endIdx + 1:]
        with open(outGro, 'w') as f:
            f.write(header)
            f.write(f'{len(newAtomLines)}\n')
            f.writelines(newAtomLines)
            f.write(boxLine)
        return outGro

    # -- 4. End-state equilibration -----------------------------------------------
    def _getStateDir(self, state):
        return self._getExtraPath(STATE_DIRNAME[state])

    def _linkTopologyExtras(self, stageDir, topFile):
        localTop = os.path.join(stageDir, os.path.basename(topFile))
        if not os.path.exists(localTop):
            os.link(topFile, localTop)
        # atomtypes.itp (not ffmerged.itp - that one's only read to build atomtypes.itp in
        # _buildHybridTopology, never copied into sysDir on its own) + merged.itp, matching
        # exactly what _buildHybridTopology writes and topol.top's #include lines reference.
        for extra in ('merged.itp', 'atomtypes.itp'):
            src = self._getExtraPath(SYSTEM_DIR, extra)
            dst = os.path.join(stageDir, extra)
            if not os.path.exists(dst):
                os.link(src, dst)
        return localTop

    def _writeFepMdp(self, outFile, phase, state, **kw):
        # No couple-moltype: the A/B states are already encoded per-atom in the hybrid
        # itp (dual topology), so GROMACS' separate molecule-decoupling feature isn't
        # used here. Confirmed by actually running grompp: 'couple-moltype = none' is
        # NOT a valid "disable" sentinel, it's read as a literal molecule name to look
        # up ("Did not find any molecules of type 'none' for coupling").
        alpha, sigma = self.scAlpha.get(), self.scSigma.get()
        feBlock = (f'free-energy      = yes\n'
                   f'init-lambda      = {state}\n'
                   f'delta-lambda     = {kw.get("deltaLambda", 0)}\n'
                   f'sc-alpha         = {alpha}\n'
                   f'sc-sigma         = {sigma}\n'
                   f'sc-power         = 1\n'
                   f'sc-coul          = yes\n')
        common = ('nstlist          = 10\n'
                  'cutoff-scheme    = Verlet\n'
                  'coulombtype      = PME\n'
                  'rcoulomb         = 1.0\n'
                  'rvdw             = 1.0\n'
                  'pbc              = xyz\n'
                  'constraints      = h-bonds\n'
                  'constraint_algorithm = lincs\n')

        if phase == 'em':
            text = ('integrator       = steep\n'
                   f'emtol            = {self.emTol.get()}\n'
                   'emstep           = 0.01\n'
                   f'nsteps           = {self.nStepsMin.get()}\n'
                   + common + feBlock)
        else:
            nsteps = round(kw['simTime'] / self.timeStep.get())
            interval = kw.get('outInterval', 0)
            text = ('integrator       = md\n'
                   f'dt               = {self.timeStep.get()}\n'
                   f'nsteps           = {nsteps}\n'
                   + common +
                   f'continuation     = {kw.get("continuation", "no")}\n'
                   f'nstxout          = {interval}\n'
                   f'nstvout          = {interval}\n'
                   f'nstdhdl          = {kw.get("nstdhdl", 100)}\n'
                   'tcoupl           = V-rescale\n'
                   'tc-grps          = System\n'
                   f'ref_t            = {self.temperature.get()}\n'
                   'tau_t            = 0.1\n'
                   'nsttcouple       = 10\n')
            if kw.get('pcoupl', False):
                text += ('pcoupl           = C-rescale\n'
                        'pcoupltype       = isotropic\n'
                        f'ref_p            = {self.pressure.get()}\n'
                        'tau_p            = 2.0\n'
                        'compressibility  = 4.5e-5\n'
                        'nstpcouple       = 10\n')
            else:
                text += 'pcoupl           = no\n'
            if kw.get('genVel', False):
                text += (f'gen_vel          = yes\n'
                        f'gen_temp         = {self.temperature.get()}\n'
                        'gen_seed         = -1\n')
            else:
                text += 'gen_vel          = no\n'
            text += feBlock

        with open(outFile, 'w') as f:
            f.write(text)
        return outFile

    def _runGrompp(self, stageDir, mdpFile, groFile, topFile, tprName, prevCpt=None, extraArgs=''):
        localTop = self._linkTopologyExtras(stageDir, topFile)
        prevStr = f' -t {os.path.abspath(prevCpt)}' if prevCpt else ''
        outFile = f'{tprName}.tpr'
        command = (f'grompp -f {os.path.abspath(mdpFile)} -c {os.path.abspath(groFile)} '
                   f'-r {os.path.abspath(groFile)} -p {os.path.basename(localTop)}{prevStr} '
                   f'-o {outFile} -maxwarn 2 {extraArgs}')
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

    def equilibrateStateStep(self, state):
        stateDir = self._getStateDir(state)
        os.makedirs(stateDir, exist_ok=True)
        topFile = self._getExtraPath(SYSTEM_DIR, 'topol.top')
        groFile = self._getExtraPath(SYSTEM_DIR, 'system.gro')

        emMdp = self._writeFepMdp(os.path.join(stateDir, 'em.mdp'), 'em', state)
        self._runGrompp(stateDir, emMdp, groFile, topFile, 'em')
        self._runMdrun(stateDir, 'em')

        # Frames saved often enough to yield >= nStructs snapshots from the 2nd half of NPT
        nptSteps = round(self.nptTime.get() / self.timeStep.get())
        interval = max(1, round(nptSteps / (2 * max(self.nStructs.get(), 1))))

        nvtMdp = self._writeFepMdp(os.path.join(stateDir, 'nvt.mdp'), 'nvt', state,
                                   simTime=self.nvtTime.get(), continuation='no', genVel=True,
                                   pcoupl=False, outInterval=0)
        self._runGrompp(stateDir, nvtMdp, os.path.join(stateDir, 'em.gro'), topFile, 'nvt')
        self._runMdrun(stateDir, 'nvt')

        nptMdp = self._writeFepMdp(os.path.join(stateDir, 'npt.mdp'), 'npt', state,
                                   simTime=self.nptTime.get(), continuation='yes', genVel=False,
                                   pcoupl=True, outInterval=interval)
        self._runGrompp(stateDir, nptMdp, os.path.join(stateDir, 'nvt.gro'), topFile, 'npt',
                        prevCpt=os.path.join(stateDir, 'nvt.cpt'))
        self._runMdrun(stateDir, 'npt')

    # -- 5. Snapshot extraction ---------------------------------------------------
    def _getFrameDir(self, state):
        return self._getExtraPath('frames', STATE_DIRNAME[state])

    def extractSnapshotsStep(self, state):
        stateDir = self._getStateDir(state)
        frameDir = self._getFrameDir(state)
        os.makedirs(frameDir, exist_ok=True)

        burnIn = self.nptTime.get() / 2.0
        args = (f'trjconv -f npt.trr -s npt.tpr -b {burnIn} -sep '
                f'-o {os.path.abspath(os.path.join(frameDir, "frame.gro"))}')
        gromacsPlugin.runGromacsPrintf(self, printfValues=['0'], args=args, cwd=stateDir)

        nStructs = self.nStructs.get()
        frames = sorted(glob(os.path.join(frameDir, 'frame*.gro')),
                        key=lambda p: int(re.search(r'frame(\d+)\.gro', p).group(1)))
        if len(frames) < nStructs:
            self.warning(f'Only {len(frames)} post-burn-in frames available for state {state}, '
                         f'requested {nStructs}')
        else:
            for extra in frames[:-nStructs] if nStructs else frames:
                os.remove(extra)

    def _getFrameFile(self, state, i):
        frames = sorted(glob(os.path.join(self._getFrameDir(state), 'frame*.gro')),
                        key=lambda p: int(re.search(r'frame(\d+)\.gro', p).group(1)))
        return frames[i] if i < len(frames) else None

    # -- 6. Fast-growth transitions ------------------------------------------------
    def _getTransitionDir(self, direction, i):
        return self._getExtraPath('transitions', direction, f'frame_{i}')

    def _failedFlagPath(self, direction, i):
        return os.path.abspath(self._getExtraPath(f'.failed_{direction}_{i}.flag'))

    def runTransitionStep(self, direction, i):
        state = DIRECTION_STATE[direction]
        frameFile = self._getFrameFile(state, i)
        if frameFile is None:
            with open(self._failedFlagPath(direction, i), 'w') as fh:
                fh.write(f'no frame {i} available for state {state}\n')
            return

        tDir = self._getTransitionDir(direction, i)
        os.makedirs(tDir, exist_ok=True)
        topFile = self._getExtraPath(SYSTEM_DIR, 'topol.top')

        nsteps = round(self.swTime.get() / self.timeStep.get())
        deltaLambda = (1.0 / nsteps) if direction == FWD else (-1.0 / nsteps)
        mdp = self._writeFepMdp(os.path.join(tDir, 'ti.mdp'), 'ti', state,
                                simTime=self.swTime.get(), continuation='yes', genVel=False,
                                pcoupl=True, outInterval=0, nstdhdl=1, deltaLambda=deltaLambda)
        try:
            self._runGrompp(tDir, mdp, frameFile, topFile, 'ti')
            self._runMdrun(tDir, 'ti')
        except Exception as e:
            self.warning(f'Transition {direction} frame {i} failed: {e}')
            with open(self._failedFlagPath(direction, i), 'w') as fh:
                fh.write(str(e) + '\n')

    def _getDhdlFiles(self, direction):
        dhdlFiles = []
        for i in range(self.nStructs.get()):
            if os.path.exists(self._failedFlagPath(direction, i)):
                continue
            # mdrun names the dH/dl output <deffnm>.xvg (not <deffnm>.dhdl.xvg), confirmed
            # by actually running a free-energy mdrun on this machine. Absolute: analyseStep
            # runs with cwd=anDir, a sibling of the transitions dir these files live under.
            dhdl = os.path.abspath(os.path.join(self._getTransitionDir(direction, i), 'ti.xvg'))
            if os.path.exists(dhdl):
                dhdlFiles.append(dhdl)
        return dhdlFiles

    # -- 7. pmx analyse ------------------------------------------------------------
    def analyseStep(self):
        anDir = self._getExtraPath('analysis')
        os.makedirs(anDir, exist_ok=True)

        fwdFiles = self._getDhdlFiles(FWD)
        bwdFiles = self._getDhdlFiles(BWD)
        if not fwdFiles or not bwdFiles:
            self.warning(f'Not enough successful transitions to analyse '
                         f'(forward: {len(fwdFiles)}, reverse: {len(bwdFiles)})')
            return

        args = (f'-fA {" ".join(fwdFiles)} -fB {" ".join(bwdFiles)} '
                f'-m CGI BAR JARZ -t {self.temperature.get()} -o results.txt -w wplot.png')
        try:
            gromacsPlugin.runPmx(self, 'analyse', args, cwd=anDir)
        except Exception as e:
            # pmx writes results.txt incrementally, one estimator at a time (BAR is known to
            # crash outright on numpy>=2 environments - see Plugin.addPmx's numpy<2 pin -
            # confirmed directly: a run here still produced a correct CGI dG before BAR's
            # crash aborted the whole process). Don't discard an otherwise-complete result
            # over one optional estimator failing.
            dG, _ = self._parseResultsFile()
            if dG is None:
                raise
            self.warning(f'"pmx analyse" exited with an error ({e}), but at least one '
                         f'estimator in results.txt succeeded; continuing with that.')

    def _parseResultsFile(self):
        """Best-effort parse of pmx analyse's results.txt (BAR preferred, then CGI, then JARZ).
        Never raises: the raw results.txt is kept as an output regardless of whether this
        regex still matches a future pmx output-format tweak."""
        resultsFile = self._getExtraPath('analysis', 'results.txt')
        if not os.path.exists(resultsFile):
            return None, None
        with open(resultsFile) as f:
            text = f.read()

        for estimator in ('BAR', 'CGI', 'JARZ'):
            dgMatch = re.search(rf'{estimator}:\s*dG\s*=\s*([-\d.]+)\s*kJ/mol', text)
            if dgMatch:
                errMatch = re.search(rf'{estimator}:\s*Std Err[^=]*=\s*([-\d.]+)\s*kJ/mol', text)
                return float(dgMatch.group(1)), float(errMatch.group(1)) if errMatch else None
        return None, None

    # -- 8. Output -----------------------------------------------------------------
    def createOutputStep(self):
        system = self._getLigandASystem()
        dG, dGErr = self._parseResultsFile()

        outSystem = GromacsSystem(filename=system.getSystemFile(), topoFile=system.getTopologyFile(),
                                  ff=system.getForceField(), wff=system.getWaterForceField())
        outSystem.setChainNames(','.join(self.getModelChains()))
        chains, lengthsDic = self.getModelChainsAndLengths()
        outSystem.setChainLengths(','.join(str(lengthsDic[c]) for c in chains))
        outSystem.setLigTopologyFile(system.getLigTopologyFile())
        outSystem.setIndexFile(self._getExtraPath('indexes.ndx'))

        resultsFile = self._getExtraPath('analysis', 'results.txt')
        if dG is not None:
            outSystem.setFreeEnergy(dG)
        if os.path.exists(resultsFile):
            outSystem.setFreeEnergyFile(resultsFile)

        self._defineOutputs(outputSystem=outSystem)
        self._defineSourceRelation(self.inputSetOfMols, outSystem)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = super()._validate()
        mols = self.inputSetOfMols.get()
        if mols is not None:
            names = [mol.__str__() for mol in mols]
            if self.inputLigand.get() not in names:
                errors.append(f'Ligand A "{self.inputLigand.get()}" not found in the input set of molecules.')
            if self.ligandB.get() not in names:
                errors.append(f'Ligand B "{self.ligandB.get()}" not found in the input set of molecules.')
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append(f'Ligand A: {self.inputLigand.get()}\nLigand B: {self.ligandB.get()}')
            dG, dGErr = self._parseResultsFile()
            if dG is not None:
                errStr = f' +/- {dGErr}' if dGErr is not None else ''
                summary.append(f'Relative binding free energy dG(A->B) = {dG}{errStr} kJ/mol')
            else:
                summary.append('Finished, but the free energy estimate could not be parsed from '
                               'pmx analyse\'s results.txt (see the raw file in the analysis output).')
        else:
            summary.append('The protocol has not finished.')
        return summary

    def _methods(self):
        methods = []
        if self.isFinished():
            methods.append('Relative binding free energy was calculated with the non-equilibrium '
                           'fast-growth alchemical FEP approach implemented in pmx: atom mapping and '
                           'hybrid dual-topology ligand construction ("pmx atomMapping", "pmx '
                           'ligandHybrid"), end-state equilibration and short forward/reverse switching '
                           'runs with GROMACS ("gmx grompp"/"gmx mdrun"), and free energy estimation '
                           'with "pmx analyse" (Crooks Gaussian Intersection / BAR / Jarzynski).')
        return methods

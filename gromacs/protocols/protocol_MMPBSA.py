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
This protocol performs MM/PBSA (or MM/GBSA) binding free energy calculations
using gmx_MMPBSA for protein–ligand complexes from GROMACS MD simulations.

Pipeline:
  1. Build a merged Protein+Ligand GROMACS index (make_ndx).
  2. Preprocess trajectory: fix PBC → fit → strip solvent/ions (trjconv × 3).
  3. Build a clean index for the stripped complex (make_ndx on com_ref.gro). ### COMPROBAR SI ES NECESARIO
  4. Write the gmx_MMPBSA input file (mmpbsa.in).
  5. Run gmx_MMPBSA.
  6. Parse results and define Scipion outputs.
"""

import os, subprocess, re, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.object import Float, String
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol

from pwchem.utils import getBaseName, convertToSdf, runOpenBabel
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem
from gromacs.protocols.protocol_MD_simulation import GromacsMDSimulation

from multiprocessing import cpu_count

from gromacs.protocols.protocol_system_prep import (
    GromacsSystemPrep,
    replaceInFile,
    GROMACS_LIST, GROMACS_MAINFF_NAME,
    GROMACS_WATERS_LIST, GROMACS_WATERFF_NAME,
    GROMACS_AMBER03, GROMACS_TIP3P,
)


# ─── Constants ────────────────────────────────────────────────────────────────
INPUT_GROMACS = 0
INPUT_MOLS    = 1

CALC_GB   = 0
CALC_PB   = 1
ENT_NONE  = 0
ENT_IE    = 1
ENT_NMODE = 2
IGB_VALS  = [1, 2, 5, 7, 8]          # AMBER igb numbers

ENTROPY_TYPE_EMPTY = "entropyType == {}"
TOPOL_TOP = "topol.top"
PREPROC_NDX = "preproc.ndx"
RESULT_DAT = "result.dat"
RESULT_CSV = "result.csv"
SCRIPT_LIGPREP_NAME = 'rdkit_addHydrogens.py'


class GromacsMmpbsa(GromacsSystemPrep):
    """
    Protein–ligand MM/PBSA / MM/GBSA binding free energy via gmx_MMPBSA.

    Mode A (Gromacs System): requires a finished MD trajectory.
    Mode B (Docked Molecules): parametrizes the ligand, builds and solvates
    the complex, minimizes it, and runs a single-frame MMPBSA calculation.

    Reference: Valdés-Tresanco et al., J. Chem. Theory Comput. 2021.
    https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/
    """
    _label = 'MM/PBSA free energy calculation'
    stepsExecutionMode = params.STEPS_PARALLEL

    def _defineParams(self, form):
        cpus = cpu_count() // 2  # don't use everything
        parallelHelp =('Threads: number of OpenMP threads assigned for minimization to fmx mdrun\n'
                        'MPIs: The number of message-passing processes assigned for calculating the binding free energy. '
                        'Usefull for trajectory calculation, not for docked molecules.\n'
                        'Total Cores Used = Scipion Threads × Max(Threads, MPIs)')
        form.addParallelSection(threads=cpus, mpi=1, binThreads=1, binThreadsHelp=parallelHelp)
        form.addSection(label='Energy calculation')
        form.addParam('inputFrom', params.EnumParam,
                      choices=['Gromacs System (post-MD)', 'Docked Molecules'],
                      default=INPUT_GROMACS,
                      label='Input type: ',
                      display=params.EnumParam.DISPLAY_HLIST,
                      help='Gromacs System: provide a finished MD run.\n'
                           'Docked Molecules: provide a docked pose set; '
                           'the protocol will prepare the system, minimize '
                           'it and run a single-frame MMPBSA calculation.')
        form.addParam('inputSetOfMols', params.PointerParam,
                        pointerClass='SetOfSmallMolecules',
                        label='Docked molecules: ',
                        condition=f'inputFrom=={INPUT_MOLS}',
                        allowsNull=True,
                        help='Set of docked molecules. The associated protein '
                             'file will be used as the receptor.')
        form.addParam('gromacsSystem', params.PointerParam,
                      label='Input Gromacs System: ',
                      pointerClass='GromacsSystem', condition=f'inputFrom=={INPUT_GROMACS}',
                      help='Gromacs system after production MD. '
                           'Must have a trajectory (.xtc), a TPR structure '
                           'file, a topology (.top) and a ligand (mol2).')

        # ── Free Energy calculation ──
        grp = form.addGroup('Calculation')
        grp.addParam('calcType', params.EnumParam,
                     label='Solvation model: ',
                     choices=['MM/GBSA', 'MM/PBSA'],
                     default=CALC_GB,
                     display=params.EnumParam.DISPLAY_COMBO,
                     help='Implicit solvation model:\n'
                          '  MM/GBSA — Generalized Born (faster, less accurate)\n'
                          '  MM/PBSA — Poisson-Boltzmann (slower, more accurate)')

        grp.addParam('igb', params.EnumParam,
                     label='GB model (igb): ',
                     choices=['1 — HCT', '2 — OBC-I', '5 — OBC-II',
                              '7 — GBn', '8 — GBn2'],
                     default=2,   # OBC-II is recommended for proteins
                     condition='calcType == {}'.format(CALC_GB),
                     help='AMBER Generalized Born model. OBC-II (igb=5) is a '
                          'good general-purpose choice for proteins.')

        grp.addParam('temperature', params.FloatParam,
                     label='Temperature (K): ',
                     default=298.15,
                     help='Specify the temperature (in K) used in the calculations.')

        grp.addParam('saltCon', params.FloatParam,
                     label='Salt concentration (M): ',
                     default=0.15,
                     help='Ionic strength for the PB calculation '
                          '(default 0.15 M, physiological). '
                          'Ions are removed from the trajectory but their '
                          'effect is modelled implicitly here.')

        # ── Entropy ──
        grp = form.addGroup('Entropy Correction')
        grp.addParam('entropyType', params.EnumParam,
                     label='Entropy method: ',
                     choices=['None',
                              'Interaction Entropy (IE)',
                              'Normal Mode Analysis (nmode)'],
                     default=ENT_NONE,
                     help='Entropy correction to obtain ΔG from ΔH:\n'
                          '  None  — skip entropy (gives ΔH; fine for ranking)\n'
                          '  IE    — fast, reasonably accurate for ranking \n'
                          '  C2    — fast, based on energy variance \n'
                          '  nmode — most rigorous but very slow')
        grp.addParam('ieSegment', params.IntParam,
                     label='IE segment (%): ',
                     default=25, expertLevel=params.LEVEL_ADVANCED,
                     condition=ENTROPY_TYPE_EMPTY.format(ENT_IE),
                     help='Percentage of frames (from the end) used to compute the '
                          'Interaction Entropy average.')
        grp = form.addGroup('Normal mode entropy calculation', condition=ENTROPY_TYPE_EMPTY.format(ENT_NMODE))

        nmodeFrame = grp.addLine('Frame selection:',  condition=ENTROPY_TYPE_EMPTY.format(ENT_NMODE),
                            help='The trajectory from which snapshots will be chosen for nmode calculations will be the collection'
                                 ' of snapshots upon which the other calculations were performed (keep low, '
                                 'e.g. 10–50 — each frame requires a minimisation). ')
        nmodeFrame.addParam('nmStartFrame', params.IntParam,
                     label='Start frame', allowsNull=True,
                     default=None,
                     help='Number of frames for normal mode entropy (keep low, '
                          'e.g. 10–50 — each frame requires a minimisation).')
        nmodeFrame.addParam('nmEndFrame', params.IntParam,
                       label='Start frame', allowsNull=True,
                       default=None,
                       help='Number of frames for normal mode entropy (keep low, '
                            'e.g. 10–50 — each frame requires a minimisation).')
        nmodeFrame.addParam('nmIntervalFrame', params.IntParam,
                       label='Interval',
                       default=1)
        nmodeMin = grp.addLine('Minimization:', condition=ENTROPY_TYPE_EMPTY.format(ENT_NMODE),
                                 help='Maximum number of minimization cycles to use per snapshot in sander. '
                                     'and convergence criteria for minimized energy gradient.')
        nmodeMin.addParam('nmMaxCycles', params.IntParam, default=10000,
                       label='Max cycles')
        nmodeMin.addParam('minConvergence', params.FloatParam, default=0.001,
                          label='Convergence')

        # ── Decomposition ──
        grp = form.addGroup('Per-Residue Decomposition')
        grp.addParam('doDecomp', params.BooleanParam,
                     label='Run decomposition analysis?: ',
                     default=False,
                     help='Decompose ΔG per residue to identify key binding '
                          'residues. Written to FINAL_DECOMP_MMPBSA.dat.')

        grp.addParam('decompResidues', params.StringParam,
                     label='Residues to print: ', default='within 6',
                     condition='doDecomp',
                     help='Which residues to include in the decomp output.\n'
                          'Examples:\n'
                          '  "within 6"   — all residues within 6 Å of the ligand (recommended)\n'
                          '  "1-10, 25"   — explicit residue range / list\n'
                          '  "all"        — every residue')

        grp.addParam('decompScheme', params.EnumParam,
                     label='Decomposition scheme: ',
                     choices=['Per residue',
                              'Per residue extended',
                              'Per atom',
                              'Per atom extended'],
                     default=0, expertLevel=params.LEVEL_ADVANCED,
                     condition='doDecomp',
                     help='1: Per-residue decomp with 1-4 terms added to internal potential terms\n'
                        '2: Per-residue decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms\n'
                        '3: Pairwise decomp with 1-4 terms added to internal potential terms\n'
                        '4: Pairwise decomp with 1-4 EEL added to EEL and 1-4 VDW added to VDW potential terms\n')

        grp.addParam('decompVerb', params.EnumParam,
                     label='Output verbosity: ', expertLevel=params.LEVEL_ADVANCED,
                     choices=['Delta (Total)', 'Delta (Decomp)', 'Full (Total)', 'Full (Decomp)'],
                     default=0, condition='doDecomp',
                     help='Level of output detail:\n'
                          '1: Delta energy totals only\n'
                          '2: Delta energy + Sidechain/Backbone breakdown\n'
                          '3: Complex/Receptor/Ligand/Delta totals\n'
                          '4: All components + Sidechain/Backbone breakdown')

        form.addSection(label='Preprocessing')
        trajPreprocGroup = form.addGroup('Trajectory preprocessing', condition=f'inputFrom=={INPUT_GROMACS}')
        self._defineTrajPreprocParams(trajPreprocGroup)

        self._defineACPYPEparams(form, condition=f'inputFrom=={INPUT_MOLS}')

        grpSys = form.addGroup('System Preparation',
                               condition=f'inputFrom=={INPUT_MOLS}')
        self._defineFFParams(grpSys)

        self._defineBoxParams(grpSys)

        grpMin = form.addGroup('Minimization', condition=f'inputFrom=={INPUT_MOLS}')
        self._defineMinimParams(grpMin)

        form.addParallelSection(threads=4, mpi=1)


    def _defineTrajPreprocParams(self, grp):
        grp.addParam('doFit', params.BooleanParam,
                     label='Fit trajectory (rot+trans)?: ',
                     default=True,
                     help='Fit trajectory to the initial structure before '
                          'stripping solvent. Strongly recommended to remove '
                          'global rotation/translation.')
        line = grp.addLine('Frame selection',
                           help='First and last frame to include in the MM/PBGB calculation. '
                                'It is common to make the calculation with the last 10ns of the simulation.')
        line.addParam('startFrame', params.IntParam,
                      label='Start: ',
                      default=1)
        line.addParam('endFrame', params.IntParam,
                      label='End: ',
                      default=1000)
        grp.addParam('interval', params.IntParam,
                     label='Frame interval: ',
                     default=10,
                     help='Take one frame every N frames (e.g. 10 = every 10th).')

    def _defineMinimParams(self, grp):
        grp.addParam('nStepsMin', params.IntParam, default=10000,
                       label='Maximum number of steps: ',
                       help='Maximum number of (minimization) steps to perform')
        grp.addParam('emTol', params.FloatParam, default=1000.0,
                       label='Maximum force objective:',
                       help='Stop minimization when the maximum force < x kJ/mol/nm.\n'
                            'https://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html#mdp-emtol')
        grp.addParam('emStep', params.FloatParam, default=0.002,
                       label='Initial step-size (nm): ',
                       help='Initial step-size (nm).\n'
                            'https://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html#mdp-emstep')

    # ── Step insertion ──────────────────────────────────────────────────────

    def _insertAllSteps(self):
        mmpbsaSteps = []  # Collect all MMPBSA calculation steps

        if self.inputFrom.get() == INPUT_MOLS:
            cStep = self._insertFunctionStep(self.convertInputStep)

            # Track parametrization steps by molName to avoid redundant parametrization
            paramSteps = {}

            for mol in self.inputSetOfMols.get():
                molName = mol.getMolName()

                # Parametrize each unique ligand once (can run in parallel for different ligands)
                if molName not in paramSteps:
                    molFile = mol.getPoseFile()
                    paramSteps[molName] = self._insertFunctionStep(
                        self.parametrizeSpecificLigandStep, molFile, molName,
                        prerequisites=[cStep])

                # Each pose follows its own pipeline
                poseId = os.path.splitext(os.path.basename(mol.getPoseFile()))[0]

                # System preparation needs its ligand to be parametrized
                sysStep = self._insertFunctionStep(
                    self.prepareSystemStep, poseId,
                    prerequisites=[paramSteps[molName]])

                # Sequential dependencies for this specific pose
                idxStep = self._insertFunctionStep(
                    self.makePreprocessingIndexStep, poseId,
                    prerequisites=[sysStep])

                preprocStep = self._insertFunctionStep(
                    self.preprocInputStep, poseId,
                    prerequisites=[idxStep])

                writeStep = self._insertFunctionStep(
                    self.writeMmpbsaInputStep, poseId,
                    prerequisites=[preprocStep])

                # MMPBSA calculation (the compute-intensive step)
                mmpbsaStep = self._insertFunctionStep(
                    self.runMmpbsaStep, poseId,
                    prerequisites=[writeStep])

                mmpbsaSteps.append(mmpbsaStep)

        else:  # INPUT_GROMACS mode - acts like serial
            idxStep = self._insertFunctionStep(self.makePreprocessingIndexStep)

            preprocStep = self._insertFunctionStep(
                self.preprocInputStep,
                prerequisites=[idxStep])

            writeStep = self._insertFunctionStep(
                self.writeMmpbsaInputStep,
                prerequisites=[preprocStep])

            mmpbsaStep = self._insertFunctionStep(
                self.runMmpbsaStep,
                prerequisites=[writeStep])

            mmpbsaSteps.append(mmpbsaStep)

        # Create output only after ALL MMPBSA calculations finish
        self._insertFunctionStep(self.createOutputStep, prerequisites=mmpbsaSteps)

    # ── Steps implementations ────────────────────────────────────────────────
    def parametrizeSpecificLigandStep(self, molFile, molName):
        """Parametrize the ligand with ACPYPE/GAFF2."""
        molFile = self.addSpecificlLigandHydrogens(molFile, molName)
        kwargs = self.getParameters()

        args = f'-i {molFile} -b {molName} -c {kwargs["chargeMethod"]} ' \
               f'-m {kwargs["multip"]} -a {kwargs["atomType"]} -q {kwargs["qprog"]} -o gmx'
        if 'netCharge' in kwargs:
            args += f' -n {kwargs["netCharge"]}'
        pwchemPlugin.runACPYPE(self, args=args, cwd=self.getLigandFileDir())


    def prepareSystemStep(self, poseId):
        """
        For one pose: pdb2gmx on receptor, merge ligand topology and
        coordinates, editconf (periodic box), solvate.
        All outputs go into the pose-specific directory.
        """
        poseDir = self.getPoseDir(poseId)
        molName = re.sub(r'_\d+$', '', poseId)
        protFile = self.getInputReceptorFile()
        sysName = os.path.splitext(os.path.basename(protFile))[0]
        mainFF = GROMACS_MAINFF_NAME[self.mainForceField.get()]
        waterFF = GROMACS_WATERFF_NAME[self.waterForceField.get()]

        # ACPYPE outputs
        ligItp = os.path.join(self.getLigParamDir(molName), f'{molName}_GMX.itp')
        ligGro = os.path.join(self.getLigParamDir(molName), f'{molName}_GMX.gro')
        # Copy ITP into pose dir so the topology #include resolves at run time
        shutil.copy(ligItp, os.path.join(poseDir, f'{molName}_GMX.itp'))

        # pdb2gmx
        print(f'Preparing mol pose: {poseId}')
        procGro = os.path.join(poseDir, f'{sysName}_processed.gro')
        p2gArgs = (f'pdb2gmx -f {protFile} -o {procGro} '
                   f'-water {waterFF} -ff {mainFF} -merge all')
        try:
            gromacsPlugin.runGromacs(self, 'gmx', p2gArgs, cwd=poseDir)
        except Exception:
            self.warning('pdb2gmx failed; retrying with -ignh')
            topFile = os.path.join(poseDir, TOPOL_TOP)
            if os.path.exists(topFile):
                os.remove(topFile)
            gromacsPlugin.runGromacs(self, 'gmx', p2gArgs + ' -ignh',
                                     cwd=poseDir)

        # merge ligand into topology and coordinates
        self.addSpecificLigandTopo(os.path.join(poseDir, TOPOL_TOP), molName)
        self.addLigandCoords(procGro, ligGro)

        # define periodic box
        boxType = self.getEnumText('boxType').lower() if self.boxType.get() != 1 else 'triclinic'
        newboxGro = os.path.join(poseDir, f'{sysName}_newbox.gro')
        ecArgs = (f'editconf -f {procGro} -o {newboxGro} '
                  f'-c -bt {boxType} ')
        ecArgs += self.getDistanceArgs()
        gromacsPlugin.runGromacs(self, 'gmx', ecArgs, cwd=poseDir)

        # solvate
        waterModel = waterFF
        if waterModel in ('spc', 'spce', 'tip3p'):
            waterModel = 'spc216'
        solvGro = os.path.join(poseDir, f'{sysName}_solv.gro')
        svArgs = (f'solvate -cp {newboxGro} -cs {waterModel}.gro '
                  f'-o {solvGro} -p {os.path.join(poseDir, TOPOL_TOP)}')
        gromacsPlugin.runGromacs(self, 'gmx', svArgs, cwd=poseDir)

        self.minimizeSystem(poseId, sysName)

    def makePreprocessingIndexStep(self, poseId=None):
        """
        Build the initial GROMACS index on the full (solvated) system.
        Creates group 21 Protein_LIG
        """
        if self.inputFrom.get() == INPUT_GROMACS:
            inputStruct = os.path.abspath(self.gromacsSystem.get().getFileName())
            cwd = self._getExtraPath()
            ndxOut = os.path.abspath(self._getExtraPath(PREPROC_NDX))
        else:
            poseDir = self.getPoseDir(poseId)
            inputStruct = os.path.abspath(os.path.join(poseDir, 'em.gro'))
            cwd = poseDir
            ndxOut = os.path.abspath(os.path.join(poseDir, PREPROC_NDX))

        args = f'make_ndx -f {inputStruct} -o {ndxOut}'
        gromacsPlugin.runGromacsPrintf(protocol=self,
                                       printfValues=['1 | 13', 'q'],
                                       args=args, cwd=cwd, mpi=False)

    def preprocInputStep(self, poseId=None):
        """
        Mode A: PBC fix + optional rot+trans fit, copy result as processTraj.xtc.
        Mode B: convert the single em.gro to a pdb
        """
        if self.inputFrom.get() == INPUT_GROMACS:
            gromacsSys = self.gromacsSystem.get()
            inputTrj   = os.path.abspath(gromacsSys.getTrajectoryFile())
            lastTpr    = os.path.abspath(gromacsSys.getTprFile())
            ndxFile    = os.path.abspath(self._getExtraPath(PREPROC_NDX))
            ligName    = gromacsSys.getLigandID()
            mergedGrp  = f'Protein_{ligName}'
            procTrj    = os.path.abspath(self._getPath('processTraj.xtc'))

            noPBC = os.path.abspath(self._getExtraPath('noPBC.xtc'))
            args  = ['trjconv', '-f', inputTrj, '-s', lastTpr,
                     '-o', noPBC, '-n', ndxFile, '-pbc', 'mol', '-center']
            gromacsPlugin.runGromacsPrintf(self,
                                           printfValues=[mergedGrp, mergedGrp],
                                           args=args, cwd=self._getExtraPath(),
                                           mpi=False)

            if self.doFit:
                fitted = os.path.abspath(self._getExtraPath('fit.xtc'))
                args   = ['trjconv', '-f', inputTrj, '-s', lastTpr,
                          '-o', fitted, '-n', ndxFile, '-fit', 'rot+trans']
                gromacsPlugin.runGromacsPrintf(self,
                                               printfValues=[mergedGrp, mergedGrp],
                                               args=args, cwd=self._getExtraPath(),
                                               mpi=False)
                shutil.copyfile(fitted, procTrj)
            else:
                shutil.copyfile(noPBC, procTrj)

        else:
            poseDir = self.getPoseDir(poseId)
            emGro   = os.path.abspath(os.path.join(poseDir, 'em.gro'))
            emPdb = os.path.abspath(os.path.join(poseDir, 'em.pdb'))

            self.convertGroToPdb(emGro, emPdb)

    def writeMmpbsaInputStep(self, poseId=None):
        """
        Generate mmpbsa.in using gmx_MMPBSA --create_input to get a fully
        annotated template, then patch the variables that the user set in the
        Scipion form via regex substitution.
        """
        inputFile = self._getMmpbsaInFile(poseId)

        # ── 1. Build --create_input keyword list ──────────────────────────
        keywords = ['gb'] if self.calcType.get() == CALC_GB else ['pb']
        if self.doDecomp.get():
            keywords.append('decomp')

        if self.entropyType.get() == ENT_NMODE:
            keywords.append('nmode')

        createArgs = ' '.join(keywords)

        gromacsPlugin.runGMXMMPBSA(self, 'gmx_MMPBSA --create_input', createArgs, cwd=os.path.dirname(inputFile))

        if not os.path.exists(inputFile):
            raise FileNotFoundError(
                f'gmx_MMPBSA --create_input did not produce mmpbsa.in file '
                f'{self._getExtraPath()}')

        with open(inputFile) as fh:
            content = fh.read()

        # ── 5. Patch variables ─────────────────────────────────
        # &general — always applied
            # Frame selection
        if self.inputFrom.get() == INPUT_MOLS:
            content = self.patchInFile(content, 'startframe', 1)
            content = self.patchInFile(content, 'endframe', 1)
            content = self.patchInFile(content, 'interval', 1)
            protName = os.path.splitext(os.path.basename(self.getInputReceptorFile()))[0]
            sysName = protName + '_' + poseId
        else:
            content = self.patchInFile(content, 'startframe', self.startFrame.get())
            content = self.patchInFile(content, 'endframe', self.endFrame.get())
            content = self.patchInFile(content, 'interval', self.interval.get())
            sysName = getBaseName(self.gromacsSystem.get().getLigTopologyFile().replace('_GMX', ''))
        content = self.patchInFile(content, 'sys_name', sysName)
        content = self.patchInFile(content, 'gmx_path', f"{gromacsPlugin.getGromacsBin(program='')}")
        content = self.patchInFile(content, 'temperature', self.temperature.get())
        content = self.patchInFile(content, 'full_traj', 1)
        content = self.patchInFile(content, 'netcdf', 1)

        # Interaction Entropy lives in &general as a flag
        if self.entropyType.get() == ENT_IE:
            content = self.patchInFile(content, 'interaction_entropy', 1)
            content = self.patchInFile(content, 'ie_segment', self.ieSegment.get())

        # &gb
        if self.calcType.get() == CALC_GB:
            igbVal = IGB_VALS[self.igb.get()]
            content = self.patchInFile(content, 'igb', igbVal)
            content = self.patchInFile(content, 'saltcon', self.saltCon.get())

        # &pb
        if self.calcType.get() == CALC_PB:
            content = self.patchInFile(content, 'istrng', f'{self.saltCon.get():.4f}')

        # decomp
        if self.doDecomp.get():
            content = self.patchInFile(content, 'idecomp',    self.decompScheme.get() + 1)
            content = self.patchInFile(content, 'dec_verbose', self.decompVerb.get())
            content = self.patchInFile(content, 'print_res',  f'"{self.decompResidues.get()}"')

        # &nmode
        if self.entropyType.get() == ENT_NMODE:
            content = self.patchInFile(content, 'nmstartframe', self.nmStartFrame.get())
            content = self.patchInFile(content, 'nmendframe', self.nmEndFrame.get())
            content = self.patchInFile(content, 'nminterval', self.nmIntervalFrame.get())
            content = self.patchInFile(content, 'maxcyc', self.nmMaxCycles.get())
            content = self.patchInFile(content, 'drms', self.minConvergence.get())

        # ── 5. Write patched file ─────────────────────────────────────────
        with open(inputFile, 'w') as fh:
            fh.write(content)

        self.info(f'mmpbsa.in ready: {inputFile}')

    def runMmpbsaStep(self, poseId=None):
        """
        Execute gmx_MMPBSA.
        Mode A: single run using the full trajectory from GromacsSystem.
        Mode B: one run per pose.
        """
        inFile  = self._getMmpbsaInFile(poseId)
        nMpi    = self.numberOfMpi.get()

        if self.inputFrom.get() == INPUT_GROMACS:
            gromacsSys   = self.gromacsSystem.get()
            lastTpr      = os.path.abspath(gromacsSys.getTprFile())
            comTrj       = os.path.abspath(self._getPath('processTraj.xtc'))
            ndxFile      = os.path.abspath(self._getExtraPath(PREPROC_NDX))
            localTopFile = os.path.abspath(self._getExtraPath('topo.top'))
            outFile      = os.path.abspath(self._getPath(RESULT_DAT))
            outCsv       = os.path.abspath(self._getPath(RESULT_CSV))
            cwd          = self._getExtraPath()

            shutil.copy(gromacsSys.getLigTopologyFile(), self._getExtraPath())
            shutil.copy(os.path.abspath(gromacsSys.getTopologyFile()), localTopFile)
        else:
            poseDir      = self.getPoseDir(poseId)
            lastTpr      = os.path.abspath(os.path.join(poseDir, 'em.tpr'))
            comTrj       = os.path.abspath(os.path.join(poseDir, 'em.pdb'))
            ndxFile      = os.path.abspath(os.path.join(poseDir, PREPROC_NDX))
            localTopFile = os.path.abspath(os.path.join(poseDir, TOPOL_TOP))
            outFile      = os.path.abspath(os.path.join(poseDir, RESULT_DAT))
            outCsv       = os.path.abspath(os.path.join(poseDir, RESULT_CSV))
            cwd          = poseDir

        args = ('-O -i {inp} -cs {cs} -ct {ct} -ci {ci} -cg 1 13 '
                '-cp {cp} -o {out} -eo {eo} -nogui').format(
            inp=inFile, cs=lastTpr, ct=comTrj,
            ci=ndxFile, cp=localTopFile, out=outFile, eo=outCsv,
        )
        gromacsPlugin.runGMXMMPBSA(self, args=args, cwd=cwd, numberOfMpi=nMpi)

    def createOutputStep(self):
        """
        Parse results and define Scipion outputs.
        Mode A: single outputSystem  (INPUT_GROMACS).
        Mode B: one annotated molecule per pose (SetOfSmallMolecules input).
        """
        calcModel = 'MMGBSA' if self.calcType.get() == CALC_GB else 'MMPBSA'

        if self.inputFrom.get() == INPUT_GROMACS:
            outFile = self._getPath(RESULT_DAT)
            if not os.path.exists(outFile):
                self.warning(f'MMPBSA output not found: {outFile}')
                return
            dg, sd = parseFinalDeltaG(outFile)
            if dg is None:
                self.warning(f'Could not parse ΔG from {outFile}')
                return
            self.info(f'{calcModel} ΔG_binding = {dg:.2f} ± {sd:.2f} kcal/mol')
            outSystem = self.gromacsSystem.get().clone()
            outSystem.setFreeEnergy(dg)
            outSystem.setFreeEnergyFile(outFile)
            self._defineOutputs(outputSystem=outSystem)

        else:
            inMols = self.inputSetOfMols.get()
            summaryLines = [f'{"Pose":<25} {"ΔG (kcal/mol)":>16}', '-' * 52]
            outputSet = inMols.createCopy(self._getPath(), copyInfo=True)

            for mol in inMols:
                annotated = self._processPoseMol(mol, calcModel, summaryLines)
                if annotated is not None:
                    outputSet.append(annotated)

            self._defineOutputs(outputSmallMolecules=outputSet)
            self.info('\n'.join(summaryLines))


    # ── Validation / info ───────────────────────────────────────────────────
    def _validate(self):
        errors = []
        if self.inputFrom.get() == INPUT_GROMACS:
            sys = self.gromacsSystem.get()
            if not sys.hasLig():
                errors.append('An input Gromacs system with Ligand is required.')
            elif not sys.getTrajectoryFile():
                errors.append('The input system must have a trajectory file (.xtc).')
        else:
            if self.numberOfMpi.get() > 1:
                errors.append('Only 1 MPI calculator is supported.')
            if self.inputSetOfMols.get() is None:
                errors.append('A set of docked molecules is required.')
        return errors

    def _summary(self):
        """Return a summary list of the final Delta G calculations."""
        calcModel = 'MM/GBSA' if self.calcType.get() == CALC_GB else 'MM/PBSA'
        summary = []

        if self.inputFrom.get() == INPUT_GROMACS:
            summary.extend(self._getGromacsSummary(calcModel))
        elif self.inputSetOfMols.get() is not None:
            summary.extend(self._getMolsSummary())

        return summary or ['Protocol has not finished yet.']

    def _getGromacsSummary(self, calcModel):
        """Helper to parse and format the summary for a GROMACS input."""
        outFile = self._getPath(RESULT_DAT)

        if not os.path.exists(outFile):
            return []

        dg, sd = parseFinalDeltaG(outFile)
        if dg is not None:
            return [f'{calcModel} ΔG_binding = {dg:.2f} ± {sd:.2f} kcal/mol']

        return []

    def _getMolsSummary(self):
        """Helper to parse and format the summary for a set of molecules."""
        summary = []

        for mol in self.inputSetOfMols.get():
            poseId = os.path.splitext(os.path.basename(mol.getPoseFile()))[0]
            outFile = os.path.join(self.getPoseDir(poseId), RESULT_DAT)

            if not os.path.exists(outFile):
                continue

            dg, _ = parseFinalDeltaG(outFile)
            if dg is not None:
                summary.append(f'Pose {poseId} ({mol.getMolName()}): {dg:.2f} kcal/mol')

        return summary

    def _methods(self):
        return [
            '{} Binding free energies were calculated with gmx_MMPBSA '
            '(Valdés-Tresanco et al., J. Chem. Theory Comput. 2021, 17, '
            '6281-6291) using the single-trajectory protocol (ST). '
            .format(
                'GBSA' if self.calcType.get() == CALC_GB else 'PBSA')
        ]

    # ── Utils ─────────────────────────────────────────────────────

    def _processPoseMol(self, mol, calcModel, summaryLines):
        """
        Process a single docked pose: parse its MMPBSA/MMGBSA result and
        annotate the molecule object.  Returns the (possibly annotated) mol
        """
        poseId = os.path.splitext(os.path.basename(mol.getPoseFile()))[0]
        outFile = os.path.join(self.getPoseDir(poseId), RESULT_DAT)

        if not os.path.exists(outFile):
            self.warning(f'Pose {poseId}: result file not found, skipping.')
            summaryLines.append(f'{poseId:<25}  {"N/A":>16}')
            return None

        dg, _ = parseFinalDeltaG(outFile)
        if calcModel == 'MMGBSA':
            mol.MMGBSA_deltaG = pwobj.Float(dg)
        else:
            mol.MMPBSA_deltaG = pwobj.Float(dg)

        summaryLines.append(f'{poseId:<25} {dg:>+16.2f}')
        return mol

    def patchInFile(self, text, key, value):
        """
        Replace the value of `key` in any namelist line, preserving
        trailing comments.  Handles:
          key                  = old_value          # comment
          key=old_value,
        """
        pattern = re.compile(
            r'^(\s*' + re.escape(key) + r'\s*=\s*)([^,#/\n]+?)(\s*(?:[,#].*)?$)',
            re.MULTILINE
        )
        newText, n = pattern.subn(r'\g<1>' + str(value) + r'\g<3>', text)
        if n == 0:
            self.warning(f'Variable "{key}" not found in template; '
                         f'it will be appended manually.')
        return newText

    def convertInputStep(self):
        sdfFiles = []
        for mol in self.inputSetOfMols.get():
            molFile = mol.getPoseFile()
            sdfFiles.append(os.path.abspath(convertToSdf(self, molFile)))

        paramFile = self.writePrepParamsFile(sdfFiles)
        pwchemPlugin.runScript(self, SCRIPT_LIGPREP_NAME, paramFile, env=RDKIT_DIC, cwd=self._getPath())

    def writePrepParamsFile(self, molFiles):
        paramsFile = self.getLigParamFile()
        with open(paramsFile, 'w') as f:
            molFilesStr = ' '.join(molFiles)
            f.write(f"ligandFiles: {molFilesStr}\n")

            f.write(f'outputDir: {self.getLigandFileDir()}\n')
            f.write('doHydrogens: True\n')
            f.write('doGasteiger: False\n')
            f.write('sanitize: False\n')
        return paramsFile

    def getLigParamFile(self):
      return os.path.abspath(self._getExtraPath('addHydrogens.txt'))

    def getLigandFileDir(self):
      lDir = os.path.abspath(self._getExtraPath('ligand'))
      if not os.path.exists(lDir):
        os.mkdir(lDir)
      return lDir

    def addSpecificlLigandHydrogens(self, inpFile, molName):
      sysbaseName = os.path.basename(inpFile).split('.')[0]

      tmpFile = os.path.abspath(self._getTmpPath(sysbaseName + '.pdb'))
      inpMol2File = os.path.abspath(self._getExtraPath(sysbaseName + '.mol2'))

      args = f'{os.path.abspath(inpFile)} -O {tmpFile}'
      runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

      args = f'{os.path.abspath(tmpFile)} -h -O {inpMol2File}'
      runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

      replaceInFile(inpMol2File, 'UNL1', 'LIG')
      return inpMol2File

    def getLigParamDir(self, molName):
        lDir = os.path.abspath(os.path.join(self.getLigandFileDir(), f'{molName}.acpype'))
        return lDir

    def getPoseDir(self, poseId):
        """Return (and create) the working directory for a given pose."""
        d = os.path.abspath(self._getExtraPath(f'poses/{poseId}'))
        os.makedirs(d, exist_ok=True)
        return d

    def addSpecificLigandTopo(self, topFile, molName):
      inStr = '\/forcefield.itp"\n'
      outStr = f'{inStr}\n; Include ligand topology\n#include "{molName}_GMX.itp"\n'
      replaceInFile(topFile, inStr, outStr)

      emptyStr = ' ' * (20-len(molName))
      inStr = '; Compound        #mols\nProtein_chain_A     1'
      outStr = f'{inStr}\n{molName}{emptyStr}1'
      replaceInFile(topFile, inStr, outStr)

    def writeIonsMDP(self, poseDir):
        outFile = os.path.join(poseDir, 'ions.mdp')
        with open(outFile, 'w') as f:
            f.write('integrator    = steep\n'
                    f'emtol         = {self.emTol.get()}\n'
                    'emstep        = 0.01\n'
                    f'nsteps        = {self.nStepsMin.get()}\n'
                    'nstlist       = 10\n'
                    'cutoff-scheme = Verlet\n'
                    'coulombtype   = cutoff\n'
                    'rcoulomb      = 1.0\n'
                    'rvdw          = 1.0\n'
                    'pbc           = xyz\n')
        return outFile

    def minimizeSystem(self, poseId, sysName):
        """Steepest-descent energy minimization for one pose."""
        poseDir = self.getPoseDir(poseId)
        topFile = os.path.join(poseDir, TOPOL_TOP)
        inputGro = self.getSolvatedGro(poseId, sysName)
        emMdp = self.writeEMMDP(poseDir)
        emTpr = os.path.join(poseDir, 'em.tpr')

        gppArgs = (f'grompp -f {emMdp} -c {inputGro} '
                   f'-p {topFile} -o {emTpr} -maxwarn 2')
        gromacsPlugin.runGromacs(self, 'gmx', gppArgs, cwd=poseDir)

        gromacsPlugin.runGromacs(self, 'gmx', f'mdrun -v -deffnm em -nt {self.binThreads.get()}',
                                 cwd=poseDir)

    def getSolvatedGro(self, poseId, sysName):
        """Path of the final solvated GRO for a pose (after ions if requested)."""
        poseDir = self.getPoseDir(poseId)
        return os.path.join(poseDir, f'{sysName}_solv.gro')

    def writeEMMDP(self, poseDir):
        outFile = os.path.join(poseDir, 'em.mdp')
        with open(outFile, 'w') as f:
            f.write('integrator    = steep\n'
                    f'emtol         = {self.emTol.get()}\n'
                    f'emstep        = {self.emStep.get()}\n'
                    f'nsteps        = {self.nStepsMin.get()}\n'
                    'nstlist       = 10\n'
                    'cutoff-scheme = Verlet\n'
                    'coulombtype   = PME\n'
                    'rcoulomb      = 1.0\n'
                    'rvdw          = 1.0\n'
                    'pbc           = xyz\n')
        return outFile

    def convertGroToPdb(self, groFile, pdbFile):
        """ Helper function to convert GRO to PDB using GROMACS editconf """
        args = f'-f {groFile} -o {pdbFile}'
        gromacsPlugin.runGromacs(self, 'gmx editconf', args)

    def _getMmpbsaInFile(self, poseId=None):
        if poseId:
            inFile = os.path.abspath(os.path.join(self.getPoseDir(poseId), 'mmpbsa.in'))
        else:
            inFile  = os.path.abspath(self._getExtraPath('mmpbsa.in'))
        return inFile

# ── helpers ─────────────────

def parseFinalDeltaG(outFile):
    """
    Extracts the final ΔTOTAL value and its Standard Deviation (SD)
    from gmx_MMPBSA output files.
    """
    pattern = re.compile(r'(?:Δ|DELTA)\s*TOTAL\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', re.IGNORECASE)

    try:
        with open(outFile, 'r', encoding='utf-8') as fh:
            for line in fh:
                m = pattern.search(line)
                if m:
                    avg = float(m.group(1))
                    sd = float(m.group(3))  # The standard SD

                    return avg, sd
    except Exception as e:
        print(f"Error reading file: {e}")

    return None, None


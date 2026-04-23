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
from pwem.protocols import EMProtocol

from pwchem.utils import getBaseName
from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem, FreeEnergyCalculation
from gromacs.protocols.protocol_MD_simulation import GromacsMDSimulation

from multiprocessing import cpu_count


# ─── Constants ────────────────────────────────────────────────────────────────
CALC_GB   = 0
CALC_PB   = 1
ENT_NONE  = 0
ENT_IE    = 1
ENT_CA    = 2
ENT_NMODE = 3
IGB_VALS  = [1, 2, 5, 7, 8]          # AMBER igb numbers shown in the enum


class GromacsMMPBSA(EMProtocol):
    """
    Protein–ligand MM/PBSA / MM/GBSA binding free energy via gmx_MMPBSA.

    The protocol prepares a consistent set of inputs (trajectory, index,
    topology, input file) and calls gmx_MMPBSA, following the recommended
    single-trajectory protocol (ST).

    Reference: Valdés-Tresanco et al., J. Chem. Theory Comput. 2021.
    https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/
    """
    _label = 'MM/PBGB binding free energy'

    def _defineParams(self, form):
        cpus = cpu_count() // 2  # don't use everything
        form.addParallelSection(mpi=1)

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('gromacsSystem', params.PointerParam,
                      label='Input Gromacs System: ',
                      pointerClass='GromacsSystem',
                      help='Gromacs system after production MD. '
                           'Must have a trajectory (.xtc), a TPR structure '
                           'file, a topology (.top) and a ligand (mol2).')

        grp = form.addGroup('Trajectory Preprocessing')
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
                     default=298.15, expertLevel=params.LEVEL_ADVANCED,
                     help='Specify the temperature (in K) used in the calculations.')

        grp.addParam('saltCon', params.FloatParam,
                     label='Salt concentration (M): ',
                     default=0.15, expertLevel=params.LEVEL_ADVANCED,
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
                              'C2 approximation',
                              'Normal Mode Analysis (nmode)'],
                     default=ENT_NONE,
                     help='Entropy correction to obtain ΔG from ΔH:\n'
                          '  None  — skip entropy (gives ΔH; fine for ranking)\n'
                          '  IE    — fast, reasonably accurate for ranking \n'
                          '  C2    — fast, based on energy variance \n'
                          '  nmode — most rigorous but very slow')

        grp = form.addGroup('Normal mode entropy calculation')
        grp.addParam('ieSegment', params.IntParam,
                     label='IE segment (%): ',
                     default=25, expertLevel=params.LEVEL_ADVANCED,
                     condition='entropyType == {}'.format(ENT_IE),
                     help='Percentage of frames (from the end) used to compute the '
                          'Interaction Entropy average.')

        nmodeFrame = grp.addLine('Frame selection:',  condition='entropyType == {}'.format(ENT_NMODE),
                            help='The trajectory from which snapshots will be chosen for nmode calculations will be the collection'
                                 ' of snapshots upon which the other calculations were performed (keep low, '
                                 'e.g. 10–50 — each frame requires a minimisation). ')

        nmodeFrame.addParam('nmStartFrame', params.IntParam,
                     label='Start frame',
                     default=1,
                     help='Number of frames for normal mode entropy (keep low, '
                          'e.g. 10–50 — each frame requires a minimisation).')
        nmodeFrame.addParam('nmEndFrame', params.IntParam,
                       label='Start frame',
                       default=10000,
                       help='Number of frames for normal mode entropy (keep low, '
                            'e.g. 10–50 — each frame requires a minimisation).')
        nmodeFrame.addParam('nmIntervalFrame', params.IntParam,
                       label='Interval',
                       default=1)
        nmodeMin = grp.addLine('Minimization:', condition='entropyType == {}'.format(ENT_NMODE),
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

    # ── Step insertion ──────────────────────────────────────────────────────
    def _insertAllSteps(self):
        self._insertFunctionStep(self.makePreprocessingIndexStep)
        self._insertFunctionStep(self.prepareTrajectoryStep)
        self._insertFunctionStep(self.writeMmpbsaInputStep)
        self._insertFunctionStep(self.runMmpbsaStep)
        self._insertFunctionStep(self.createOutputStep)

    # ── Step implementations ────────────────────────────────────────────────
    def makePreprocessingIndexStep(self):
        """
        Build the initial GROMACS index on the full (solvated) system.
        Creates group 21 Protein_LIG
        """
        inputStruct = os.path.abspath(self.gromacsSystem.get().getFileName())
        ndxOut = self._getExtraPath('preproc.ndx')  # keep for reference elsewhere

        args = "make_ndx -f {} -o {}".format(inputStruct, os.path.abspath(ndxOut))

        gromacsPlugin.runGromacsPrintf(protocol=self, printfValues=["1 | 13", "q"], args=args,
            cwd=self._getExtraPath(), mpi=False)

    def prepareTrajectoryStep(self):
        """
        Trajectory preprocessing pipeline:
          trjconv 1 — fix PBC (-pbc mol -center), output System
          trjconv 2 — fit to initial structure (-fit rot+trans)  [optional]
          trjconv 3 — strip to Protein+Ligand only
          trjconv 4 — extract reference frame 0 as com_ref.gro
        """
        gromacsSys    = self.gromacsSystem.get()
        inputTrj   = os.path.abspath(gromacsSys.getTrajectoryFile())
        inputStruct = os.path.abspath(gromacsSys.getFileName())
        lastTpr = os.path.abspath(gromacsSys.getTprFile())
        ndxFile    = os.path.abspath(self._getExtraPath('preproc.ndx'))
        ligName    =gromacsSys.getLigandID()
        mergedGrp  = 'Protein_{}'.format(ligName)

        # Step 1: Fix PBC (center on receptor, output all atoms)
        noPBC = os.path.abspath(self._getExtraPath('noPBC.xtc'))
        args = [
            "trjconv",
            "-f", inputTrj,
            "-s", lastTpr,  # Notice we are pointing to the full .tpr here!
            "-o", noPBC,
            "-n", ndxFile,
            "-pbc", "mol",
            "-center"
        ]

        cmd2 = [mergedGrp, mergedGrp]
        gromacsPlugin.runGromacsPrintf(self, printfValues=[mergedGrp, mergedGrp], args=args,
                                       cwd=self._getExtraPath(), mpi=False)

        #  3 (optional): Fit rot+trans
        if self.doFit:
            fitted = os.path.abspath(self._getExtraPath('fit.xtc'))
            args = [
                "trjconv",
                "-f", inputTrj,
                "-s", lastTpr,
                "-o", fitted,
                "-n", ndxFile,
                "-fit", "rot+trans"
            ]
            cmd2 = [mergedGrp, mergedGrp]

            gromacsPlugin.runGromacsPrintf(self, printfValues=[mergedGrp, mergedGrp], args=args,
                                           cwd=self._getExtraPath(), mpi=False)
            preFinalTrj = fitted
        else:
            preFinalTrj = noPBC

        procTrj = os.path.abspath(self._getPath('processTraj.xtc'))
        shutil.copyfile(preFinalTrj, procTrj)

    def makeMmpbsaIndexStep(self):
        """
        Build the index file for the *stripped* complex (com_ref.gro).
        In the stripped system only Protein and Ligand atoms remain, so the
        default make_ndx groups already contain both. We just save it and
        verify the expected groups are present.
        """
        comRef  = self._getExtraPath('com_ref.gro')
        ndxOut  = self._getExtraPath('mmpbsa.ndx')
        recGrp  = self.receptorGroup.get().strip()
        ligName = self.ligandName.get().strip()

        gromacsPlugin.runGromacsPrintf(
            printfValues=['q'],
            args=' make_ndx -f {} -o {}'.format(comRef, ndxOut),
            cwd=self._getExtraPath()
        )

        # Persist IDs so runMmpbsaStep can use them without re-parsing
        self._storeGroupIds(recIdx, ligIdx)
        self.info('MMPBSA index: receptor group {} (id={}), '
                  'ligand group {} (id={})'.format(recGrp, recIdx, ligName, ligIdx))

    def writeMmpbsaInputStep(self):
        """
        Generate mmpbsa.in using gmx_MMPBSA --create_input to get a fully
        annotated template, then patch the variables that the user set in the
        Scipion form via regex substitution.
        """
        inputFile = self._getExtraPath('mmpbsa.in')

        # ── 1. Build --create_input keyword list ──────────────────────────
        keywords = ['gb'] if self.calcType.get() == CALC_GB else ['pb']
        if self.doDecomp.get():
            keywords.append('decomp')

        if self.entropyType.get() == ENT_NMODE:
            keywords.append('nmode')
        # IE entropy is controlled via a &general variable, not a separate

        createArgs = ' '.join(keywords)

        # print(f'Generating input template: gmx_MMPBSA --create_input {createArgs}')
        gromacsPlugin.runGMXMMPBSA(self, 'gmx_MMPBSA --create_input', createArgs, cwd=self._getExtraPath())

        if not os.path.exists(inputFile):
            raise FileNotFoundError(
                f'gmx_MMPBSA --create_input did not produce mmpbsa.in file '
                f'{self._getExtraPath()}')

        with open(inputFile) as fh:
            content = fh.read()

        # ── 5. Patch variables with regex ─────────────────────────────────
        # &general — always applied
        content = self.patchInFile(content, 'startframe', self.startFrame.get())
        content = self.patchInFile(content, 'endframe', self.endFrame.get())
        content = self.patchInFile(content, 'interval', self.interval.get())
        content = self.patchInFile(content, 'sys_name', getBaseName(self.gromacsSystem.get().getLigTopologyFile().replace("_GMX", "")))
        content = self.patchInFile(content, 'gmx_path', f"{gromacsPlugin.getGromacsBin(program='')}")
        content = self.patchInFile(content, 'temperature', self.temperature.get())
        content = self.patchInFile(content, 'full_traj', 1)
        content = self.patchInFile(content, 'netcdf', 1)

        # Interaction Entropy lives in &general as a flag
        if self.entropyType.get() == ENT_IE:
            content = self.patchInFile(content, 'interaction_entropy', 1)
            content = self.patchInFile(content, 'ie_segment', self.ieSegment.get())
        elif self.entropyType.get() == ENT_CA:
            content = self.patchInFile(content, 'c2_entropy', 1)

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

    def runMmpbsaStep(self):
        """
        Execute gmx_MMPBSA with the prepared inputs.
        gmx_MMPBSA must be available in the active environment (PATH).
        """
        extra   = self._getExtraPath
        gromacsSys = self.gromacsSystem.get()
        comTrj  = os.path.abspath(self._getPath('processTraj.xtc'))
        ndxFile = os.path.abspath(extra('preproc.ndx'))
        inFile  = os.path.abspath(extra('mmpbsa.in'))
        lastTpr = os.path.abspath(gromacsSys.getTprFile())

        topFile = os.path.abspath(self.gromacsSystem.get().getTopologyFile())
        localTopFile = os.path.abspath(self._getExtraPath('topo.top'))
        outFile = os.path.abspath(self._getPath('result.dat'))
        outCsv = os.path.abspath(self._getPath('result.csv'))

        shutil.copy(gromacsSys.getLigTopologyFile(), self._getExtraPath())
        shutil.copy(topFile, localTopFile)

        args = ('-O '
                '-i {inp} '
                '-cs {cs} '
                '-ct {ct} '
                '-ci {ci} '
                '-cg 1 13 '
                '-cp {cp} '
                '-o {out} '
                '-eo {outCsv} '
                '-nogui').format(
            inp=inFile, cs=lastTpr, ct=comTrj,
            ci=ndxFile, cp=localTopFile, out=outFile, outCsv=outCsv
        )
        # Build the same command string but force bash
        condaHook = '/home/joaquin/miniconda/etc/profile.d/conda.sh'
        envName = 'gmxMMPBSA-1.6.4'
        nMpi = self.numberOfMpi.get()
        mmpbsaCmd = f'mpirun -np {nMpi} gmx_MMPBSA {args}'

        fullCmd = (
            f'bash -c "source {condaHook} && '
            f'conda activate {envName} && '
            f'{mmpbsaCmd}"'
        )

        gromacsPlugin.runGMXMMPBSA(self, args=args, cwd=self._getExtraPath(), numberOfMpi=nMpi)

    def createOutputStep(self):
        # ...
        outFile = self._getPath('result.dat')
        outCsv = self._getPath('result.csv')

        if not os.path.exists(outFile):
            self.warning('MMPBSA output file not found: {}'.format(outFile))
            return

        dg, sd = _parseFinalDeltaG(outFile)

        if dg is not None:
            calcModel = 'MM/GBSA' if self.calcType.get() == CALC_GB else 'MM/PBSA'

            # Update the print statement to include SD
            self.infoAG = '{} ΔG_binding = {:.2f} ± {:.2f} kcal/mol'.format(calcModel, dg, sd)
            print(self.infoAG)

            output = FreeEnergyCalculation(
                deltaG=dg,
                deviation=sd,
                calcType=calcModel,
                resultFile=outFile,
                csvFile=outCsv,
                entropyType=self.entropyType.get()
            )

            self._defineOutputs(outputFreeEnergy=output)

        else:
            self.warning('Could not parse ΔG_binding from {}'.format(outFile))

    # ── Validation / info ───────────────────────────────────────────────────
    def _validate(self):
        errors = []
        gromacsSys = self.gromacsSystem.get()
        if gromacsSys is None:
            errors.append('An input Gromacs system is required.')
            return errors
        if not gromacsSys.getTrajectoryFile():
            errors.append('The input system must have a trajectory file (.xtc).')
        return errors

    def _summary(self):
        summary = [f'Results in {self.getPath("result.dat")}']
        outFile = self._getPath('result.dat')
        dg, sd = _parseFinalDeltaG(outFile)
        model = 'MM/GBSA' if self.calcType.get() == CALC_GB else 'MM/PBSA'
        summary.append(f'{model} ΔG_binding = {dg:.2f} ± {sd:.2f} kcal/mol"')
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

    def _loadGroupIds(self):
        """Load receptor/ligand group IDs written by makeMmpbsaIndexStep."""
        idFile = self._getExtraPath('_group_ids.txt')
        with open(idFile) as fh:
            lines = fh.read().strip().split()
        return int(lines[0]), int(lines[1])

    def runInteractiveCommand(self, program, args, inputs):
        """
        Pipes interactive inputs into a program and executes via Scipion's runJob.

        :param program: The binary to run (e.g., 'gmx')
        :param args: The command arguments as a string (e.g., 'make_ndx -f sys.gro')
        :param inputs: List of strings to send (e.g., ['1 | 13', 'q'])
        """
        # Join inputs with \n for the printf command
        inputString = "\\n".join(str(i) for i in inputs) + "\\n"

        fullProgram = f'printf "{inputString}" | {program}'
        self.runJob(fullProgram, args)

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
        new_text, n = pattern.subn(r'\g<1>' + str(value) + r'\g<3>', text)
        if n == 0:
            self.warning(f'Variable "{key}" not found in template; '
                         f'it will be appended manually.')
        return new_text


# ── helpers ─────────────────

def _parseFinalDeltaG(outFile):
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
                    sd_prop = float(m.group(2))  # Optional, if you need it
                    sd = float(m.group(3))  # The standard SD

                    return avg, sd
    except Exception as e:
        print(f"Error reading file: {e}")

    return None, None

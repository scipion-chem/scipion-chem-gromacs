# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


"""
This module will prepare the system for the simumlation
"""
import os
import shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem


class GromacsModifySystem(EMProtocol):
    """
    AI Generated:

    Protocol to modify the structure and trajectory of a GROMACS system.

    It applies the usual trajectory post-processing operations with gmx trjconv and
    gmx filter, producing a new GromacsSystem ready for analysis or visualization.

    Overview
    --------
    The available operations, applied in this order, are:

    1. PBC correction of the trajectory
    2. Cleaning: keep only the protein (and ligand), removing waters and ions
    3. Fitting to the reference structure
    4. Cutting between two times
    5. Subsampling of the frames
    6. Low-pass or high-pass filtering

    Input
    -----
    gromacsSystem:
        GromacsSystem with the structure and (optionally) the trajectory to modify.

    correctPBC:
        Three-step trjconv correction: whole molecules, continuous path (nojump), and
        centering on the protein (or protein-ligand complex) with compact wrapping.
        Requires the .tpr file of the system.

    cleaning:
        Keep only the protein, or the protein-ligand complex when the system has a
        ligand (Protein_<ligand> group of the system index).

    doFit / fitting:
        Fit the frames to the reference structure (rot+trans, rotxy+transxy,
        translation, transxy or progressive).

    doDrop / firstTime, lastTime, timeUnit:
        Keep only the frames between two times (0 means from the first / until the last).

    doSubsample / subsampleF:
        Keep one frame every N.

    doFiltering / filter, filterF:
        Low-pass (with its filter length) or high-pass filtering of the trajectory.

    Workflow
    --------
    1. modifySystem:
       - Extracts the cleaned structure with editconf, if requested
       - Corrects the PBC jumps, if requested
       - Runs a single trjconv with the cleaning, fitting, cutting and subsampling options,
         using the index file of the system
       - Runs gmx filter, if requested

    2. createOutputStep:
       - Creates the new GromacsSystem and reads the time and frame range of its trajectory

    Output
    ------
    outputSystem:
        GromacsSystem with the modified structure and trajectory and the topology of the
        input system.

    Key Features
    ------------
    - Standard PBC correction procedure
    - Protein-ligand aware cleaning
    - Fitting, cutting, subsampling and filtering in one protocol

    Notes
    -----
    - The low-pass filter length must be at least 2
    - If the input system has no trajectory, only the structure is modified
    """
    _label = 'System modification'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('gromacsSystem', params.PointerParam, label="Input Gromacs System: ",
                      pointerClass='GromacsSystem',
                      help='Gromacs solvated system to be simulated')

        group = form.addGroup('PBC correction')
        group.addParam('correctPBC', params.BooleanParam, label="Correct PBC jumps?: ", default=False,
                       help='Corrects periodic boundary condition (PBC) artifacts in the trajectory before '
                            'any other processing (cleaning, fitting, cutting, subsampling, filtering). '
                            'Runs the standard 3-step gmx trjconv procedure:\n\n'
                            '1) "-pbc whole": re-assembles molecules split across opposite sides of the '
                            'box within a frame (fixes intra-molecular bonds).\n'
                            '2) "-pbc nojump": builds a smooth, continuous path through space over time, '
                            'using the whole molecules from step 1 (fixes inter-frame teleportation).\n'
                            '3) "-center -pbc mol -ur compact": centers the protein and neatly wraps solvent/ions '
                            'around it. \n\n')

        group = form.addGroup('Cleaning')
        group.addParam('cleaning', params.BooleanParam, label="Clean protein?: ", default=False,
                       help='Remove waters and ions from the system, keeping only the protein')
        group = form.addGroup('Fitting')
        group.addParam('doFit', params.BooleanParam, label="Fit trajectory?: ", default=False,
                       help='Fit trajectory to initial structure')
        group.addParam('fitting', params.EnumParam, label="Fitting type: ", default=0, condition='doFit',
                       choices=['rot+trans', 'rotxy+transxy', 'translation', 'transxy', 'progressive'],
                       help='Fitting technique to the initial structure. '
                            'https://manual.gromacs.org/documentation/5.1/onlinehelp/gmx-trjconv.html')
        group = form.addGroup('Cutting')
        group.addParam('doDrop', params.BooleanParam, label="Cut trajectory?: ", default=False,
                       help='Cut a trajectory, saving only from first to last times')
        line = group.addLine('Times: ', condition='doDrop',
                             help='First and last times to save from trajectory (0 from first, 0 until last)')
        line.addParam('firstTime', params.FloatParam, label="First: ", default=0)
        line.addParam('lastTime', params.FloatParam, label="Last: ", default=0)
        line.addParam('timeUnit', params.EnumParam, label="Units: ",
                      default=1, choices=['fs', 'ps', 'ns', 'us', 'ms', 's'])

        group = form.addGroup('Subsample')
        group.addParam('doSubsample', params.BooleanParam, label="Subsample trajectory?: ", default=False,
                       help='Subsample trajectory frames')
        group.addParam('subsampleF', params.IntParam, label="Subsample factor: ", default=10, condition='doSubsample',
                       help='Subsample factor. Take a frame for each x original frames')

        group = form.addGroup('Filtering')
        group.addParam('doFiltering', params.BooleanParam, label="Filter trajectory?: ", default=False,
                       help='Perform filter on trajectory frames')
        group.addParam('filter', params.EnumParam, label="Filter trajectory: ", default=0,
                       choices=['Low pass', 'High pass'], condition='doFiltering',
                       help='Type of filter to use in the trajectory, specially useful if you subsample')
        group.addParam('filterF', params.IntParam, label="Filter length: ", default=10,
                       condition='doFiltering and filter==0',
                       help='Sets the filter length as well as the output interval for low-pass filtering')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.modifySystem)
        self._insertFunctionStep(self.createOutputStep)

    def modifySystem(self):
        system = self.gromacsSystem.get()
        inputStructure = os.path.abspath(system.getFileName())
        inputTrajectory = system.getTrajectoryFile()
        indexFile = gromacsPlugin.ensureIndexFile(self)
        outputGroup = system.getComplexGroup() if self.cleaning.get() else 'System'

        if self.cleaning:
            params = " editconf -f {} -n {} -o {}".format(inputStructure, indexFile, self.getCleanStructureFile())
            gromacsPlugin.runGromacsPrintf(self, printfValues=[outputGroup],
                                           args=params, cwd=self._getPath())
            if system.getTprFile():
                # A .tpr with only the kept atoms, as the one of the input has all the waters and ions
                params = " convert-tpr -s {} -n {} -o {}".format(os.path.abspath(system.getTprFile()), indexFile,
                                                                self.getCleanTprFile())
                gromacsPlugin.runGromacsPrintf(self, printfValues=[outputGroup],
                                               args=params, cwd=self._getPath())
        else:
            shutil.copy(inputStructure, self.getCleanStructureFile())

        if inputTrajectory:
            inputTrajectory = os.path.abspath(inputTrajectory)
            if self.correctPBC:
                inputTrajectory = self.correctPBCJumps(inputTrajectory)

            auxTrj = os.path.abspath(self._getExtraPath('cleanTrajectory.xtc'))
            convArgs = " trjconv -f {} -s {} -n {} -o {}". \
                format(inputTrajectory, inputStructure, indexFile, auxTrj)

            extraArgs = ''
            if self.doFit:
                extraArgs += ' -fit {}'.format(self.getEnumText('fitting'))
            if self.doDrop:
                firstTime, lastTime = self.getCutTime(self.firstTime.get()), self.getCutTime(self.lastTime.get())
                if firstTime != 0:
                    extraArgs += ' -b {}'.format(firstTime)
                if lastTime != 0:
                    extraArgs += ' -e {}'.format(lastTime)
                if firstTime or lastTime:
                    extraArgs += ' -tu {}'.format(self.getEnumText('timeUnit'))

            if self.doSubsample:
                extraArgs += ' -skip {}'.format(self.subsampleF.get())

            if extraArgs or self.cleaning:
                convArgs += extraArgs
                # trjconv asks for the fitting group (only with -fit) and then for the output group
                printfVals = ['Protein', outputGroup] if self.doFit else [outputGroup]
                gromacsPlugin.runGromacsPrintf(self, printfValues=printfVals,
                                               args=convArgs, cwd=self._getPath())
            else:
                auxTrj = inputTrajectory

            if self.doFiltering:
                filterArgs = 'filter -f {}'.format(auxTrj)
                if self.filter.get() == 0:
                    filterArgs += ' -nf {}'.format(self.filterF.get())
                filterArgs += self.getFilteringArgs(self.getCleanTrajectoryFile(), self.getCleanStructureFile())
                if '-fit' in filterArgs:
                    gromacsPlugin.runGromacsPrintf(self, printfValues=['Protein'],
                                                   args=filterArgs, cwd=self._getPath())
                else:
                    gromacsPlugin.runGromacs(self, args=filterArgs, cwd=self._getPath())
            else:
                shutil.copy(auxTrj, self.getCleanTrajectoryFile())

    def createOutputStep(self):
      inSystem = self.gromacsSystem.get()
      outSystem = GromacsSystem()
      outSystem.setSystemFile(os.path.relpath(self.getCleanStructureFile()))
      outSystem.setTopologyFile(inSystem.getTopologyFile())
      if inSystem.getTrajectoryFile():
          outSystem.setTrajectoryFile(os.path.relpath(self.getCleanTrajectoryFile()))
          outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())

      if inSystem.getMinimizedFile():
          outSystem.setMinimizedFile(inSystem.getMinimizedFile())
      if inSystem.getTprFile():
          outSystem.setTprFile(os.path.relpath(self.getCleanTprFile()) if self.cleaning else inSystem.getTprFile())
      outSystem.setLigandID(inSystem.getLigandID())
      outSystem.setLigTopologyFile(inSystem.getLigTopologyFile())
      outSystem.setChainNames(inSystem.getChainNames())
      outSystem.setChainLengths(inSystem.getChainLengths())

      self._defineOutputs(outputSystem=outSystem)


    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
        if self.doFiltering and self.filter.get() == 0 and self.filterF.get() <= 1:
            vals.append('The factor for low-pass filtering needs to be at least 2')

        if self.correctPBC:
            tprFile = self.gromacsSystem.get().getTprFile()
            if not tprFile or not os.path.exists(os.path.abspath(tprFile)):
                vals.append('A .tpr file is necessary to correct PBC jumps.')

        return vals

    def _warnings(self):
        warns = []
        inputTrajectory = self.gromacsSystem.get().getTrajectoryFile()
        if not inputTrajectory:
            if self.correctPBC:
                warns.append('The input system has no trajectory, so no PBC correction will be performed')
            if self.doFit:
                warns.append('The input system has no trajectory, so no fitting will be performed')
            if self.doDrop:
                warns.append('The input system has no trajectory, so no cutting will be performed')
            if self.doSubsample:
                warns.append('The input system has no trajectory, so no subsampling will be performed')
            if self.doFiltering:
                warns.append('The input system has no trajectory, so no filtering will be performed')

        return warns

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def getCleanStructureFile(self):
        inputStructure = self.gromacsSystem.get().getFileName()
        name, ext = os.path.splitext(inputStructure)
        return os.path.abspath(self._getPath(os.path.basename(name)) + ext)

    def getCleanTprFile(self):
        return os.path.abspath(self._getPath('cleanSystem.tpr'))

    def getCleanTrajectoryFile(self):
        inputTrajectory = self.gromacsSystem.get().getTrajectoryFile()
        return os.path.abspath(self._getPath(os.path.basename(inputTrajectory.split(".")[0])) + '.xtc')

    def getFilteringArgs(self, trjFile, strFile=''):
        if self.getEnumText('filter') == 'Low pass':
            args = ' -ol {}'.format(trjFile)
        elif self.getEnumText('filter') == 'High pass':
            args = ' -oh {} -s {}'.format(trjFile, strFile)
            if not self.doFit:
                #Trajectory needs to be fitted when using high pass filtering. If not done previously
                args += ' -fit'
                print('Do fitting for high pass filtering')
        return args

    def getCutTime(self, time):
        intTime = 0 if (time.is_integer() and int(time) == 0) else time
        return intTime

    def correctPBCJumps(self, inputTrajectory):
        """Fixes periodic boundary condition (PBC) artifacts of the raw trajectory following the
        strict 3-step gmx trjconv procedure: whole molecules -> continuous (nojump) path ->
        centered / compactly wrapped system. The centering group is the protein-ligand complex
        if the system has a parametrized ligand, otherwise the protein.
        """
        inpSystem = self.gromacsSystem.get()
        indexFile = gromacsPlugin.ensureIndexFile(self)
        centerGroup = inpSystem.getComplexGroup()

        # -pbc mol requires a .tpr for -s
        tprFile = inpSystem.getTprFile()
        tprFile = os.path.abspath(tprFile) if tprFile else None
        if not tprFile or not os.path.exists(tprFile):
            raise FileNotFoundError('A .tpr file is necessary to correct PBC jumps, but the input Gromacs '
                                    'System has none associated.')

        wholeTrj = os.path.abspath(self._getTmpPath('pbcWhole.xtc'))
        wholeArgs = ' trjconv -f {} -s {} -n {} -pbc whole -o {}'. \
            format(inputTrajectory, tprFile, indexFile, wholeTrj)
        gromacsPlugin.runGromacsPrintf(self, printfValues=['System'], args=wholeArgs, cwd=self._getPath())

        nojumpTrj = os.path.abspath(self._getTmpPath('pbcNojump.xtc'))
        nojumpArgs = ' trjconv -f {} -s {} -n {} -pbc nojump -o {}'. \
            format(wholeTrj, tprFile, indexFile, nojumpTrj)
        gromacsPlugin.runGromacsPrintf(self, printfValues=['System'], args=nojumpArgs, cwd=self._getPath())

        correctedTrj = os.path.abspath(self._getExtraPath('pbcCorrected.xtc'))
        centerArgs = ' trjconv -f {} -s {} -n {} -center -pbc mol -ur compact -o {}'. \
            format(nojumpTrj, tprFile, indexFile, correctedTrj)
        gromacsPlugin.runGromacsPrintf(self, printfValues=[centerGroup, 'System'],
                                       args=centerArgs, cwd=self._getPath())

        return correctedTrj
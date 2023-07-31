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
    This protocol modifies a gromacs system trajectory and/or coordinates:
        - Cleans from waters and ions
        - Subsamples trajectory and applies filters
        - Fits trajectory to initial structure
    """
    _label = 'system modification'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('gromacsSystem', params.PointerParam, label="Input Gromacs System: ",
                      pointerClass='GromacsSystem',
                      help='Gromacs solvated system to be simulated')
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
        self._insertFunctionStep('modifySystem')
        self._insertFunctionStep('createOutputStep')

    def modifySystem(self):
        inputStructure = os.path.abspath(self.gromacsSystem.get().getFileName())
        inputTrajectory = self.gromacsSystem.get().getTrajectoryFile()

        if self.cleaning:
            params = " make_ndx -f {} -o clean.ndx".format(os.path.abspath(inputStructure))
            gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'q'],
                                           args=params, cwd=self._getPath())

            params = " editconf -f {} -n clean.ndx -o {}".format(inputStructure, self.getCleanStructureFile())
            gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'q'],
                                           args=params, cwd=self._getPath())
        else:
            shutil.copy(inputStructure, self.getCleanStructureFile())

        if inputTrajectory:
            inputTrajectory = os.path.abspath(inputTrajectory)
            auxTrj = os.path.abspath(self._getExtraPath('cleanTrajectory.xtc'))
            convArgs = " trjconv -f {} -s {} -o {}". \
                format(inputTrajectory, inputStructure, auxTrj)
            
            extraArgs = ''
            if self.cleaning:
                extraArgs += ' -n clean.ndx'
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

            if extraArgs:
                convArgs += extraArgs
                gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'Protein'],
                                               args=convArgs, cwd=self._getPath())
            else:
                auxTrj = inputTrajectory

            if self.doFiltering:
                filterArgs = 'filter -f {}'.format(auxTrj)
                if self.filter.get() == 0:
                    filterArgs += ' -nf {}'.format(self.filterF.get())

                filterArgs += self.getFilteringArgs(self.getCleanTrajectoryFile(), self.getCleanStructureFile())

                if '-fit' in filterArgs:
                    gromacsPlugin.runGromacsPrintf(printfValues=['Protein'],
                                                   args=filterArgs, cwd=self._getPath())
                else:
                    gromacsPlugin.runGromacs(self, args=filterArgs, cwd=self._getPath())
            else:
                shutil.copy(auxTrj, self.getCleanTrajectoryFile())



    def createOutputStep(self):
      outSystem = GromacsSystem()
      outSystem.setOriStructFile(self.getCleanStructureFile())
      outSystem.setSystemFile(self.getCleanStructureFile())
      outSystem.setTopologyFile(self.gromacsSystem.get().getTopologyFile())
      if self.gromacsSystem.get().getTrajectoryFile():
          outSystem.setTrajectoryFile(os.path.relpath(self.getCleanTrajectoryFile()))
          outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())

      self._defineOutputs(outputSystem=outSystem)


    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
        if self.doFiltering and self.filter.get() == 0 and self.filterF.get() <= 1:
            vals.append('The factor for low-pass filtering needs to be at least 2')

        return vals

    def _warnings(self):
        warns = []
        inputTrajectory = self.gromacsSystem.get().getTrajectoryFile()
        if not inputTrajectory:
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
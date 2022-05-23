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

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem


class GromacsCleanSystem(EMProtocol):
    """
    This protocol cleans a gromacs system trajectory and/or coordinates from waters and ions
    """
    _label = 'system cleaning'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('gromacsSystem', params.PointerParam, label="Input Gromacs System: ",
                      pointerClass='GromacsSystem',
                      help='Gromacs solvated system to be simulated')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('cleanSystem')
        self._insertFunctionStep('createOutputStep')

    def cleanSystem(self):
        inputStructure = self.gromacsSystem.get().getFileName()
        inputTrajectory = self.gromacsSystem.get().getTrajectoryFile()

        params_grompp = " make_ndx -f {} -o clean.ndx".format(os.path.abspath(inputStructure))
        gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'q'],
                                       args=params_grompp, cwd=self._getPath())

        params_grompp = " editconf -f {} -n clean.ndx -o {}".\
          format(os.path.abspath(inputStructure), self.getCleanStructureFile())
        gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'q'],
                                       args=params_grompp, cwd=self._getPath())

        if inputTrajectory:
            params_grompp = " trjconv -f {} -n clean.ndx -o {}".\
              format(os.path.abspath(inputTrajectory), self.getCleanTrajectoryFile())
            gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'q'],
                                           args=params_grompp, cwd=self._getPath())

    def createOutputStep(self):
      outSystem = GromacsSystem(filename=self.getCleanStructureFile())
      outSystem.setTopologyFile(self.gromacsSystem.get().getTopologyFile())
      if self.gromacsSystem.get().getTrajectoryFile():
          outSystem.setTrajectoryFile(self.getCleanTrajectoryFile())

      self._defineOutputs(outputSystem=outSystem)


    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
        return vals

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def getCleanStructureFile(self):
        inputStructure = self.gromacsSystem.get().getFileName()
        return os.path.abspath(self._getPath(os.path.basename(inputStructure.split(".")[0])) + '.gro')

    def getCleanTrajectoryFile(self):
        inputTrajectory = self.gromacsSystem.get().getTrajectoryFile()
        return os.path.abspath(self._getPath(os.path.basename(inputTrajectory.split(".")[0])) + '.xtc')


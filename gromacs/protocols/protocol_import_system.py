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
This module will import a system from the topology, coordinates [and trajectory] files
"""
import os
import shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from gromacs import Plugin as gromacsPlugin
from gromacs.objects import GromacsSystem


class GromacsImportSystem(EMProtocol):
    """
    This protocol import a gromacs system trajectory and/or coordinates:
    """
    _label = 'import system'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputCoords', params.FileParam, label="Input Gromacs Coordinates (gro): ", allowsNull=False,
                      help='Gromacs coordinates file (gro)')
        form.addParam('inputTopology', params.FileParam, label="Input Gromacs Topology (top): ", allowsNull=False,
                      help='Gromacs topology file (top)')
        form.addParam('inputTrajectory', params.FileParam, label="Input Gromacs Trajectory (trr, xtc): ",
                      help='Gromacs trajectory file (xtc / trr)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
      outSystem = GromacsSystem()
      outSystem.setSystemFile(self.inputCoords.get())
      outSystem.setOriStructFile(self.inputCoords.get())
      outSystem.setTopologyFile(self.inputTopology.get())
      if self.inputTrajectory.get():
          outSystem.setTrajectoryFile(self.inputTrajectory.get())
          outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())

      self._defineOutputs(outputSystem=outSystem)


    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
        return vals

    def _warnings(self):
        warns = []
        return warns

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

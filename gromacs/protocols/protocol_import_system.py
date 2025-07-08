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

User IA Manual: ImportSystem Protocol

The ImportSystem protocol is used to bring into Scipion-Chem a molecular system
that has been previously prepared for simulation with GROMACS. This system
typically includes a set of coordinate and topology files, which together define
the atomic structure, molecular interactions, and simulation box configuration.
By importing the system into the workflow, users can integrate it with downstream
protocols for energy minimization, molecular dynamics, or analysis.

To run the protocol, the user must provide at least the structure file in GRO
format, which contains the atomic coordinates and box vectors. In addition, a
topology file is required to define the molecular components, parameters, and
force field assignments. This file is usually in TOP format and may reference
other files such as ITP or include directives. The protocol ensures that all
referenced files are accessible and properly parsed for integration.

The user may also specify whether the imported system contains solvent, ions, or
restraints, which can influence how the system is treated in later steps. For
example, solvent molecules may be needed for pressure coupling, while positional
restraints may be required to maintain the stability of certain components during
equilibration. The protocol records all such metadata to ensure reproducibility
and correct behavior in downstream simulations.

Once executed, the protocol registers the molecular system as a Scipion object,
preserving the directory structure and file dependencies. This object becomes
the reference for subsequent GROMACS protocols, including energy minimization,
equilibration, and production runs. The imported system is not modified, and
users retain full control over how it is used in the workflow.

In summary, the ImportSystem protocol serves as the entry point for GROMACS-ready
systems in Scipion-Chem. It enables seamless integration of externally prepared
systems into structured workflows for molecular simulation, ensuring consistency,
traceability, and compatibility with other components of the platform.

    """
    _label = 'import system'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputCoords', params.FileParam, label="Input Gromacs Coordinates (gro, pdb): ", allowsNull=False,
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

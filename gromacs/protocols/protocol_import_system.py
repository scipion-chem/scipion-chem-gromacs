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

    AI Generated:

        GromacsImportSystem

        Overview
        --------
        This protocol imports a molecular system prepared for GROMACS simulations
        into Scipion-Chem.

        It registers coordinate, topology, and optional trajectory files as a
        structured system object that can be used in downstream molecular dynamics
        workflows.

        The protocol acts as the entry point for GROMACS-based simulations within
        Scipion.

        Inputs
        ------
        inputCoords:
            Coordinate file defining atomic positions and simulation box.
            Supported formats:
            - GRO
            - PDB

        inputTopology:
            Topology file describing molecular structure, force field parameters,
            and system composition.
            Typically in TOP format (may include ITP dependencies).

        inputTrajectory:
            Optional trajectory file containing simulation frames.
            Supported formats:
            - XTC
            - TRR

        Workflow
        --------
        1. Input acquisition
           - Reads coordinate and topology files
           - Optionally reads trajectory file

        2. System construction
           - Initializes a GromacsSystem object
           - Assigns:
             - coordinate file
             - topology file
             - original structure reference

        3. Trajectory integration (optional)
           - If trajectory is provided:
             - Registers trajectory file
             - Extracts trajectory metadata (frames, timing, etc.)

        4. Output registration
           - Stores system as a Scipion-compatible object
           - Preserves file references and structure

        Output
        ------
        outputSystem:
            GromacsSystem object containing:
            - Coordinates (structure)
            - Topology (force field and connectivity)
            - Optional trajectory data
            - Simulation metadata

        Summary
        -------
        This protocol imports a pre-built GROMACS system into Scipion,
        enabling:
        - integration with molecular dynamics workflows
        - reproducible simulation setup
        - reuse of externally prepared systems
        - compatibility with energy minimization and MD protocols

        Notes
        -----
        - Does not modify input files (read-only import)
        - Requires consistent topology and coordinate files
        - External topology dependencies (e.g., ITP files) must be accessible
        - Trajectory input is optional but recommended for analysis workflows

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
        self._insertFunctionStep(self.createOutputStep)

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

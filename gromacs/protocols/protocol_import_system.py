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
import re
import shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import AtomicStructHandler

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
        form.addParam('inputCoords', params.FileParam, label="Input Gromacs Coordinates (gro, pdb): ", allowsNull=False,
                      help='Gromacs coordinates file (gro or pdb)')
        form.addParam('inputTopology', params.FileParam, label="Input Gromacs Topology (top): ", allowsNull=False,
                      help='Gromacs topology file (top)')
        form.addParam('inputTopologyIncludeDir', params.PathParam,
                      label="Topology include directory: ", allowsNull=True, default='',
                      help='Optional directory with files included inside the topology file, '
                           'by #include.')
        form.addParam('inputTrajectory', params.FileParam, label="Input Gromacs Trajectory (trr, xtc): ",
                      help='Gromacs trajectory file (xtc / trr)')
        form.addParam('inputIndex', params.FileParam, label="Input Gromacs Index (ndx): ",
                      allowsNull=True,
                      help='Optional Gromacs index file (ndx). It will be copied and used by downstream '
                           'simulation protocols when creating or selecting index groups.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
      outSystem = GromacsSystem()
      coordsFile = self.inputCoords.get()
      localCoords = self.copyToProtocolDir(coordsFile)
      chainNames = ','.join(self.getModelChains())

      if not chainNames:
          chainNames = ','.join(self.inferChainNamesGro(localCoords))

      systemCoords = self.convertPdbToGro(localCoords)
      outSystem.setSystemFile(systemCoords)
      outSystem.setOriStructFile(localCoords)

      topologyFile = self.inputTopology.get()
      localTop = self.copyToProtocolDir(topologyFile)
      includeDir = self.inputTopologyIncludeDir.get()
      if includeDir:
          localIncludeDir = self._getPath(os.path.basename(includeDir.rstrip(os.sep)))
          if os.path.abspath(includeDir) != os.path.abspath(localIncludeDir):
              shutil.copytree(includeDir, localIncludeDir, dirs_exist_ok=True)
          outSystem.setTopologyIncludeDir(localIncludeDir)
      outSystem.setTopologyFile(localTop)

      if self.inputTrajectory.get():
          localTrajectory = self.copyToProtocolDir(self.inputTrajectory.get())
          outSystem.setTrajectoryFile(localTrajectory)
          outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())
      if self.inputIndex.get():
          gromacsPlugin.createIndex
          localIndex = self.copyToProtocolDir(self.inputIndex.get())
          outSystem.setIndexFile(localIndex)

      if chainNames:
          outSystem.setChainNames(chainNames)

      self._defineOutputs(outputSystem=outSystem)

    def copyToProtocolDir(self, inputFile):
        localFile = self._getPath(os.path.basename(inputFile))
        if os.path.abspath(inputFile) != os.path.abspath(localFile):
            shutil.copyfile(inputFile, localFile)
        return localFile

    def inferChainNamesGro(self, groFile):
        chains = []
        with open(groFile) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')) and len(line) > 21:
                    chain = line[21].strip()
                    if chain and chain not in chains:
                        chains.append(chain)
        return chains

    def convertPdbToGro(self, coordsFile):
        if os.path.splitext(coordsFile)[1].lower() not in ['.pdb', '.ent']:
            return coordsFile

        groFile = self._getPath(os.path.splitext(os.path.basename(coordsFile))[0] + '.gro')
        if not os.path.exists(groFile):
            args = '-f {} -o {}'.format(os.path.abspath(coordsFile), os.path.abspath(groFile))
            gromacsPlugin.runGromacs(self, 'gmx editconf', args, cwd=self._getPath())
        return groFile

    def getModelChains(self):
        inputStructure = self.inputCoords.get()
        if not inputStructure.endswith('.pdb'):
          return self.inferChainNamesGro(inputStructure)

        structureHandler = AtomicStructHandler()
        structureHandler.read(inputStructure)
        structureHandler.getStructure()
        chains, _ = structureHandler.getModelsChains()
        return list(chains[0].keys())

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
        includeDir = self.inputTopologyIncludeDir.get()
        if includeDir and not os.path.isdir(includeDir):
            vals.append('Topology include directory does not exist or is not a directory: {}'.format(includeDir))
        elif includeDir:
            includeBase = os.path.basename(includeDir.rstrip(os.sep))
            for includePath in self.getQuotedTopologyIncludes():
                normPath = os.path.normpath(includePath)
                if os.path.isabs(normPath):
                    continue
                parts = normPath.split(os.sep)
                if parts[0] == includeBase:
                    copiedPath = os.path.join(includeDir, *parts[1:])
                    if not os.path.exists(copiedPath):
                        vals.append('Topology include not found in selected directory: {}'.format(includePath))
        if self.inputIndex.get() and not os.path.isfile(self.inputIndex.get()):
            vals.append('Input index file does not exist or is not a file: {}'.format(self.inputIndex.get()))
        return vals

    def getQuotedTopologyIncludes(self):
        includes = []
        includeRe = re.compile(r'^\s*#\s*include\s+"([^"]+)"')
        with open(self.inputTopology.get()) as f:
            for line in f:
                match = includeRe.match(line)
                if match:
                    includes.append(match.group(1))
        return includes

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

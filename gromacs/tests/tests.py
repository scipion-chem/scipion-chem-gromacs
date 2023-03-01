# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from gromacs.protocols import GromacsSystemPrep, GromacsModifySystem, GromacsMDSimulation
from gromacs import Plugin as gromacsPlugin

workflow = '''{'simTime': 100.0, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': False, 'trajInterval': 1.0, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'Protein', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'Energy min', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman'}
{'simTime': 0.1, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': True, 'trajInterval': 0.05, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'MainChain', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman'}
{'simTime': 0.2, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': True, 'trajInterval': 0.05, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'None', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'NPT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman'}
'''
summary = '''1) Minimization (steep): 100 steps, 1000.0 objective force, restraint on Protein, 300.0 K
2) MD simulation: 0.1 ps, NVT ensemble, restraint on MainChain, 300.0 K
3) MD simulation: 0.2 ps, NPT ensemble, 300.0 K'''



class TestGromacsPrepareSystem(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPrepareSystem(cls):
        protPrepare = cls.newProtocol(
            GromacsSystemPrep,
            inputStructure=cls.protImportPDB.outputPdb,
            boxType=1, sizeType=1, padDist=2.0,
            mainForceField=0, waterForceField=2,
            placeIons=1, cationType=7, anionType=1)

        cls.launchProtocol(protPrepare)
        return protPrepare

    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


class TestGromacsRunSimulation(TestGromacsPrepareSystem):

    def _runSimulation(self, protPrepare):
        protSim = self.newProtocol(
            GromacsMDSimulation,
            gromacsSystem=protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)

        outFile = self.createGroupsFile(protPrepare.outputSystem, inIndex=None, outIndex=protSim.getCustomIndexFile(),
                                        outFile=protSim.getCustomGroupsFile())
        self.launchProtocol(protSim)
        return protSim

    def createGroupsFile(self, system, inIndex=None, outIndex='/tmp/indexes.ndx', outFile='/tmp/indexGroups.txt',
                           inputCommands=['q']):
        outDir = os.path.dirname(outFile)
        inIndex = ' -n {}'.format(inIndex) if inIndex else ''
        command = 'make_ndx -f {}{} -o {} > {}'.format(system.getSystemFile(), inIndex, outIndex, outFile)

        if not inputCommands[-1] == 'q':
            inputCommands.append('q')
        gromacsPlugin.runGromacsPrintf(printfValues=inputCommands, args=command, cwd=outDir)
        return outFile

    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSim, 'outputSystem', None))


class TestGromacsTrajMod(TestGromacsRunSimulation):

    def _modSimulation(self, protSim):
        protMod = self.newProtocol(
            GromacsModifySystem,
            gromacsSystem=protSim.outputSystem, cleaning=True, doFit=True)

        self.launchProtocol(protMod)
        return protMod

    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        protMod = self._modSimulation(protSim)
        self.assertIsNotNone(getattr(protMod, 'outputSystem', None))

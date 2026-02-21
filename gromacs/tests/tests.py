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
import numpy as np

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb, exists

from gromacs.protocols import (GromacsSystemPrep, GromacsModifySystem, GromacsMDSimulation,
                               PlumedRunAnalysis)
from gromacs.protocols.protocol_plumed_alone import ALPHARMSD, DISTANCE, ANGLE, TORSION
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
        protSim.setObjLabel('gromacs - gmx MD sim')

        outIndex = protSim.getCustomIndexFile()
        if os.path.exists(outIndex):
            protSim.parseIndexFile(outIndex)
        else:
            protSim.createIndexFile(protPrepare.outputSystem, inIndex=None, outIndex=protSim.getCustomIndexFile())

        self.launchProtocol(protSim)
        return protSim
    
    def _runSimulationMPI(self, protPrepare):
        protSim = self.newProtocol(
            GromacsMDSimulation, gmxMPI=True, numberOfMpi=2, 
            gromacsSystem=protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)
        protSim.setObjLabel('gromacs - gmx_mpi MD sim')

        outIndex = protSim.getCustomIndexFile()
        if os.path.exists(outIndex):
            protSim.parseIndexFile(outIndex)
        else:
            protSim.createIndexFile(protPrepare.outputSystem, inIndex=None, outIndex=protSim.getCustomIndexFile())

        self.launchProtocol(protSim)
        return protSim

    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSim, 'outputSystem', None))
        protSimMPI = self._runSimulationMPI(protPrepare)
        self._waitOutput(protSimMPI, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSimMPI, 'outputSystem', None))


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


class TestPlumedAnalysis(BaseTest):
    @classmethod
    def setUpClass(cls):
        # Create a new project
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
        cls._runImportPdb2()
        cls.protPrepare = cls._runPrepareSystem()
        cls._waitOutput(cls.protPrepare, 'outputSystem', sleepTime=10)
        cls.protSim = cls._runSimulation(cls.protPrepare)
        cls._waitOutput(cls.protSim, 'outputSystem', sleepTime=10)

    def testPlumedDistance(self):
        protPlumed1 = self.newProtocol(PlumedRunAnalysis, measureType=DISTANCE)
        protPlumed1.inputSystem.set(self.protImportPdb3o21.outputPdb)
        protPlumed1.setObjLabel('Plumed_distance_A')
        self.launchProtocol(protPlumed1)

        self._waitOutput(protPlumed1, 'outputSystem', sleepTime=100)
        self.assertIsNotNone(getattr(protPlumed1, 'outputSystem', None))

        self.assertTrue(exists(protPlumed1._getPath('COLVAR')))

        with open(protPlumed1._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 1), 3.0)


    def testPlumedAngle(self):
        protPlumed2 = self.newProtocol(PlumedRunAnalysis, measureType=ANGLE)
        protPlumed2.inputSystem.set(self.protImportPdb3o21.outputPdb)
        protPlumed2.selection1.set('chain A and name CA and resid 46')
        protPlumed2.selection2.set('chain A and name CA and resid 116')
        protPlumed2.selection3.set('chain A and name CA and resid 196')
        protPlumed2.setObjLabel('Plumed_angle_A')
        self.launchProtocol(protPlumed2)

        self._waitOutput(protPlumed2, 'outputSystem', sleepTime=100)
        self.assertIsNotNone(getattr(protPlumed2, 'outputSystem', None))

        self.assertTrue(exists(protPlumed2._getPath('COLVAR')))

        with open(protPlumed2._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 2), 1.24)


    def testPlumedTorsion(self):
        protPlumed3 = self.newProtocol(PlumedRunAnalysis, measureType=TORSION)
        protPlumed3.inputSystem.set(self.protImportPdb3o21.outputPdb)
        protPlumed3.setObjLabel('Plumed_torsion_AB')
        self.launchProtocol(protPlumed3)

        self._waitOutput(protPlumed3, 'outputSystem', sleepTime=100)
        self.assertIsNotNone(getattr(protPlumed3, 'outputSystem', None))

        self.assertTrue(exists(protPlumed3._getPath('COLVAR')))

        with open(protPlumed3._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 2), -0.57)


    def testPlumedTorsion2(self):
        protPlumed4 = self.newProtocol(PlumedRunAnalysis, measureType=TORSION)
        protPlumed4.inputSystem.set(self.protImportPdb3o21.outputPdb)
        protPlumed4.setObjLabel('Plumed_torsion_CD')
        protPlumed4.selection1.set('chain C and name CA and resid 117 to 243,chain C and name CA and resid 354 to 380')
        protPlumed4.selection2.set('chain C and name CA and resid 4 to 116,chain C and name CA and resid 244 to 353')
        protPlumed4.selection3.set('chain D and name CA and resid 4 to 116,chain D and name CA and resid 244 to 353')
        protPlumed4.selection4.set('chain D and name CA and resid 117 to 243,chain D and name CA and resid 354 to 380')
        self.launchProtocol(protPlumed4)

        self._waitOutput(protPlumed4, 'outputSystem', sleepTime=100)
        self.assertIsNotNone(getattr(protPlumed4, 'outputSystem', None))

        self.assertTrue(exists(protPlumed4._getPath('COLVAR')))

        with open(protPlumed4._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 2), -0.34)


    def testPlumedOnPreppedSystem(self):
        protPlumed5 = self.newProtocol(PlumedRunAnalysis, measureType=DISTANCE)
        protPlumed5.inputSystem.set(self.protPrepare.outputSystem)
        protPlumed5.setObjLabel('Plumed_AK_prep_dist')
        protPlumed5.selection1.set('chain A and name CA and resid 14')
        protPlumed5.selection2.set('chain A and name CA and resid 131')
        self.launchProtocol(protPlumed5)
        self.assertIsNotNone(getattr(protPlumed5, 'outputSystem', None))

        self.assertTrue(exists(protPlumed5._getPath('COLVAR')))

        with open(protPlumed5._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 1), 0.8)
        self.assertEqual(len(lines[1:]), 1)

    def testPlumedOnRunSystem(self):
        protPlumed5 = self.newProtocol(PlumedRunAnalysis, measureType=DISTANCE)
        protPlumed5.inputSystem.set(self.protSim.outputSystem)
        protPlumed5.setObjLabel('Plumed_AK_sim_dist')
        protPlumed5.selection1.set('chain A and name CA and resid 14')
        protPlumed5.selection2.set('chain A and name CA and resid 131')
        self.launchProtocol(protPlumed5)
        self.assertIsNotNone(getattr(protPlumed5, 'outputSystem', None))

        self.assertTrue(exists(protPlumed5._getPath('COLVAR')))

        with open(protPlumed5._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[1].strip().split()[-1]), 1), 0.8)
        self.assertEqual(len(lines[1:]), protPlumed5.outputSystem.getNumFrames())


    def testPlumedAlpharmsdHelix(self):
        protPlumed1 = self.newProtocol(PlumedRunAnalysis, measureType=ALPHARMSD)
        protPlumed1.inputSystem.set(self.protImportPDB.outputPdb)
        protPlumed1.setObjLabel('Plumed_alpharmsd_helix')
        protPlumed1.selection1.set('resid 65-70')
        self.launchProtocol(protPlumed1)

        self._waitOutput(protPlumed1, 'outputSystem', sleepTime=100)
        self.assertIsNotNone(getattr(protPlumed1, 'outputSystem', None))

        self.assertTrue(exists(protPlumed1._getPath('COLVAR')))

        with open(protPlumed1._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 1), 1.0)

    def testPlumedAlpharmsdStrand(self):
        protPlumed1 = self.newProtocol(PlumedRunAnalysis, measureType=ALPHARMSD)
        protPlumed1.inputSystem.set(self.protImportPDB.outputPdb)
        protPlumed1.setObjLabel('Plumed_alpharmsd_strand')
        protPlumed1.selection1.set('resid 105-110')
        self.launchProtocol(protPlumed1)

        self._waitOutput(protPlumed1, 'outputSystem', sleepTime=100)
        self.assertIsNotNone(getattr(protPlumed1, 'outputSystem', None))

        self.assertTrue(exists(protPlumed1._getPath('COLVAR')))

        with open(protPlumed1._getPath('COLVAR'), 'r') as fi:
            lines = fi.readlines()

        self.assertEqual(np.round(float(lines[-1].strip().split()[-1]), 1), 0.0)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_mut1.pdb'))
        cls.protImportPDB.setObjLabel('Input PDB 1ake mut1')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runImportPdb2(cls):
        cls.protImportPdb3o21 = cls.newProtocol(ProtImportPdb, inputPdbData=0,
                                        pdbId='3o21')
        cls.protImportPdb3o21.setObjLabel('Input PDB 3o21')
        cls.launchProtocol(cls.protImportPdb3o21)

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

    @classmethod
    def _runSimulation(cls, protPrepare):
        protSim = cls.newProtocol(
            GromacsMDSimulation,
            gromacsSystem=protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)
        protSim.setObjLabel('gromacs - gmx MD sim')

        outIndex = protSim.getCustomIndexFile()
        if os.path.exists(outIndex):
            protSim.parseIndexFile(outIndex)
        else:
            protSim.createIndexFile(protPrepare.outputSystem, inIndex=None, outIndex=protSim.getCustomIndexFile())

        cls.launchProtocol(protSim)
        return protSim

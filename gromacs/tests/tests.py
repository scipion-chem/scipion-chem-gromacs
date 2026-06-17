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

import os, logging, tempfile

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from gromacs.protocols import GromacsSystemPrep, GromacsModifySystem, GromacsMDSimulation, GromacsMmpbsa
from gromacs import Plugin as gromacsPlugin

from pwchem.tests import TestExtractLigand
from pwchem.protocols.VirtualDrugScreening.protocol_receptor_preparation import ProtChemPrepareReceptor

STRUCTURE, LIGAND = 0, 1
chainStr = '{"model": 0, "chain": "A", "residues": 236}'

workflow = '''{'simTime': 100.0, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': False, 'trajInterval': 1.0, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'Protein', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'Energy min', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman', 'coupleStyle': 'isotropic'}
{'simTime': 0.1, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': True, 'trajInterval': 0.05, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'MainChain', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman', 'coupleStyle': 'isotropic'}
{'simTime': 0.2, 'timeStep': 0.002, 'nStepsMin': 100, 'emStep': 0.002, 'emTol': 1000.0, 'timeNeigh': 10, 'saveTrj': True, 'trajInterval': 0.05, 'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1, 'restraints': 'None', 'restraintForce': 50.0, 'integrator': 'steep', 'ensemType': 'NPT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman', 'coupleStyle': 'isotropic'}
'''
summary = '''1) Minimization (steep): 100 steps, 1000.0 objective force, restraint on Protein, 300.0 K
2) MD simulation: 0.1 ps, NVT ensemble, restraint on MainChain, 300.0 K
3) MD simulation: 0.2 ps, NPT ensemble, 300.0 K'''


class TestGromacsPrepareSystem(TestExtractLigand):

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='1uaz')
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareReceptor(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            HETATM=True, rchains=True,
            chain_name=chainStr)

        cls.launchProtocol(cls.protPrepareReceptor)

    @classmethod
    def _runPrepareSystem(cls, protPrepare, inputFrom=STRUCTURE):
        protPrepareS = cls.newProtocol(
            GromacsSystemPrep, inputFrom=inputFrom)

        if inputFrom == STRUCTURE:
            protPrepareS.inputStructure.set(protPrepare)
            protPrepareS.inputStructure.setExtended('outputStructure')
        else:
            protPrepareS.inputSetOfMols.set(protPrepare)
            protPrepareS.inputSetOfMols.setExtended('outputSmallMolecules')
            protPrepareS.inputLigand.set('SmallMolecule (g1_1uaz_RET-1_1 molecule)')

        cls.launchProtocol(protPrepareS)
        return protPrepareS

    def test(self):
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)

        protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))

    def test2(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


class TestGromacsRunSimulation(TestGromacsPrepareSystem):

    def _runSimulation(self, protPrepare):
        protSim = self.newProtocol(
            GromacsMDSimulation,
            gromacsSystem=protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)
        protSim.setObjLabel('gromacs - gmx MD sim')

        self.launchProtocol(protSim)
        return protSim

    def _runSimulationMPI(self, protPrepare):
        protSim = self.newProtocol(
            GromacsMDSimulation, gmxMPI=True, numberOfMpi=2,
            gromacsSystem=protPrepare.outputSystem, workFlowSteps=workflow, summarySteps=summary)
        protSim.setObjLabel('gromacs - gmx_mpi MD sim')

        self.launchProtocol(protSim)
        return protSim

    def test(self):
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)

        protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSim, 'outputSystem', None))
        protSimMPI = self._runSimulationMPI(protPrepare)
        self._waitOutput(protSimMPI, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protSimMPI, 'outputSystem', None))

    test2 = None


class TestGromacsTrajMod(TestGromacsRunSimulation):

    def _modSimulation(self, protSim):
        protMod = self.newProtocol(
            GromacsModifySystem,
            gromacsSystem=protSim.outputSystem, cleaning=True, doFit=True)

        self.launchProtocol(protMod)
        return protMod

    def test(self):
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)

        protPrepare = self._runPrepareSystem(self.protPrepareReceptor)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)
        protMod = self._modSimulation(protSim)
        self.assertIsNotNone(getattr(protMod, 'outputSystem', None))
    test2 = None


class TestGromacsMMPBSA(TestGromacsRunSimulation, TestExtractLigand):
    @classmethod
    def _runInteractions(cls, protIn, inputFrom=STRUCTURE):
        protInt = cls.newProtocol(
          GromacsMmpbsa, inputFrom=inputFrom, nStepsMin=500, interval=1)

        if inputFrom == STRUCTURE:
            protInt.gromacsSystem.set(protIn)
            protInt.gromacsSystem.setExtended('outputSystem')
        else:
            protInt.inputSetOfMols.set(protIn)
            protInt.inputSetOfMols.setExtended('outputSmallMolecules')

        cls.launchProtocol(protInt)
        return protInt

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protPrepare = self._runPrepareSystem(protExtract, inputFrom=LIGAND)
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)

        protSim = self._runSimulation(protPrepare)
        self._waitOutput(protSim, 'outputSystem', sleepTime=10)

        protInt = self._runInteractions(protSim)
        self._waitOutput(protInt, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protInt, 'outputSystem', None))

    def test2(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protInt = self._runInteractions(protExtract, inputFrom=LIGAND)
        self._waitOutput(protInt, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protInt, 'outputSmallMolecules', None))


class TestGromacsRenumberResidues(BaseTest):
    """Unit test for the continuous residue renumbering performed during system preparation.
    It does not require the GROMACS binaries, only Biopython, since it exercises
    GromacsSystemPrep.renumberResiduesPDB directly."""

    @staticmethod
    def _writeTwoChainPDB(path):
        """Two chains (A, B) that BOTH restart their residue numbering at 1, reproducing the situation
        where merging the chains yields a .gro with repeated residue numbers."""
        def atom(serial, name, resname, chain, resseq, x, y, z, elem):
            return (f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain}{resseq:>4}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {elem:>2}\n")

        lines, serial = [], 1
        for chain in ('A', 'B'):
            for resseq, resname in [(1, 'PRO'), (2, 'GLN'), (3, 'ILE')]:
                for name, elem in [('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O')]:
                    lines.append(atom(serial, name, resname, chain, resseq,
                                       serial * 0.1, serial * 0.1, serial * 0.1, elem))
                    serial += 1
            lines.append("TER\n")
        lines.append("END\n")
        with open(path, 'w') as f:
            f.writelines(lines)

    @staticmethod
    def _chainResNums(pdb):
        from Bio import PDB
        structure = PDB.PDBParser(QUIET=True).get_structure('x', pdb)
        return {ch.id: [r.id[1] for r in ch] for ch in list(structure)[0]}

    def test(self):
        from gromacs.protocols.protocol_system_prep import GromacsSystemPrep

        inPdb = os.path.join(tempfile.mkdtemp(), 'two_chain_dup.pdb')
        outPdb = inPdb.replace('_dup.pdb', '_renum.pdb')
        self._writeTwoChainPDB(inPdb)

        # before: both chains numbered 1-3 (duplicates across chains)
        before = self._chainResNums(inPdb)
        self.assertEqual(before['A'], [1, 2, 3])
        self.assertEqual(before['B'], [1, 2, 3])

        # run the real method with a minimal object (it only needs self._log)
        class _Fake:
            _log = logging.getLogger('test')
        GromacsSystemPrep.renumberResiduesPDB(_Fake(), inPdb, outPdb)

        after = self._chainResNums(outPdb)
        allNums = [n for nums in after.values() for n in nums]
        # no duplicates, continuous from the original start, chains preserved
        self.assertEqual(len(allNums), len(set(allNums)), 'Duplicate residue numbers remain')
        self.assertEqual(allNums, list(range(1, len(allNums) + 1)))
        self.assertEqual(after['A'], [1, 2, 3])
        self.assertEqual(after['B'], [4, 5, 6])

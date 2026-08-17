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

from gromacs.protocols import GromacsSystemPrep, GromacsModifySystem, GromacsMDSimulation, GromacsMmpbsa, \
    GromacsPmxRBFE, GromacsPmxABFE
from gromacs import Plugin as gromacsPlugin

from pwchem.tests import TestExtractLigand
from pwchem.protocols.VirtualDrugScreening.protocol_receptor_preparation import ProtChemPrepareReceptor

STRUCTURE, LIGAND = 0, 1
chainStr = '{"model": 0, "chain": "A", "residues": 236}'
# JNK1 chain A (data/tests/smallMolecules/FEP/jnk1_18625-1_18626-1.pdb), residues 1-358 -
# real residue range confirmed against the source PDB used to build that file.
jnk1ChainStr = '{"model": 0, "chain": "A", "residues": 358}'

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
            protPrepareS.inputLigand.set('SmallMolecule (g1_1uaz_RET_255-1_1 molecule)')

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


class TestGromacsPmxRBFE(TestGromacsPrepareSystem):
    """Smoke test for the pmx relative binding free energy (RBFE) protocol.

    `test` uses a real published RBFE benchmark edge instead of docking anything:
    ligands 18625-1/18626-1 against JNK1, from pmx's own protLig_benchmark
    (https://github.com/deGrootLab/pmx/tree/master/protLig_benchmark) and its
    ligand_tutorial.ipynb, both grounded in Gapsys/Perez-Benito et al. 2020,
    Chem. Sci. 11:1140 (the paper pmx's non-equilibrium RBFE approach is
    validated against). Deliberately NOT built via docking (this plugin's tests
    must not depend on a separate docking plugin like autodock): the two
    ligands' real, already co-located poses (confirmed: ligand centroid ~2 A
    from a real JNK1 chain-A atom) were taken directly from the tutorial's own
    input files and merged with JNK1 chain A into one PDB
    (data/tests/smallMolecules/FEP/jnk1_18625-1_18626-1.pdb), extracted here via
    the same ProtExtractLigands mechanism the self-transform test below already
    uses - no docking step anywhere in this test. The two ligands are also
    genuinely similar (RDKit Morgan Tanimoto 0.77 - one chlorine moved to a
    different ring position, exactly the kind of edge RBFE is meant for),
    unlike an earlier version of this test that blind-docked unrelated ZINC
    compounds (Tanimoto <0.26) into 1uaz.

    Experimental ddG for this edge is -3.22 kJ/mol (from the tutorial) - not
    asserted here since the settings below are tiny/fast for testing, not for
    a physically meaningful free energy estimate, but available as a real
    target if this test is ever run with production-scale settings.

    `test2` keeps the original self->self (RET->RET) transformation on 1uaz,
    which should give dG(A->B)~=0 within error - a sanity check on the
    alchemical machinery itself, independent of ligand-pair identity.

    See claude/decisions/pmx_RBFE.md for the full comparison against the
    tutorial/paper's methodology (in particular: this protocol currently only
    computes the bound-complex leg, not the solvent/water leg the published
    method subtracts to get a true binding free energy - flagged there, not
    fixed here)."""

    @classmethod
    def _runImportJNK1PDB(cls):
        dsLig = DataSet.getDataSet('smallMolecules')
        cls.protImportJNK1PDB = cls.newProtocol(
            ProtImportPdb, inputPdbData=1,
            pdbFile=dsLig.getFile('FEP/jnk1_18625-1_18626-1.pdb'))
        cls.launchProtocol(cls.protImportJNK1PDB)
        return cls.protImportJNK1PDB

    @classmethod
    def _runPmxRBFE(cls, protExtract):
        # getMolName() is guessed from the shared source PDB's filename and comes back
        # identical for both ligands extracted from the same file - the pose/file path
        # (which keeps ProtExtractLigands' own per-residue naming, "..._L25_900.cif" vs
        # "..._L26_901.cif") is what actually distinguishes them here.
        molA = str(next(m for m in protExtract.outputSmallMolecules if 'L25' in m.getPoseFile()))
        molB = str(next(m for m in protExtract.outputSmallMolecules if 'L26' in m.getPoseFile()))

        protRBFE = cls.newProtocol(
            GromacsPmxRBFE, nStructs=2, swTime=2.0, nvtTime=2.0, nptTime=4.0, nStepsMin=200)
        protRBFE.inputSetOfMols.set(protExtract)
        protRBFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protRBFE.inputLigand.set(molA)
        protRBFE.ligandB.set(molB)
        protRBFE.setObjLabel('gromacs - pmx RBFE (JNK1 18625-1 -> 18626-1)')

        cls.launchProtocol(protRBFE)
        return protRBFE

    @classmethod
    def _runPmxRBFESelfTransform(cls, protExtract):
        protRBFE = cls.newProtocol(
            GromacsPmxRBFE, nStructs=2, swTime=2.0, nvtTime=2.0, nptTime=4.0, nStepsMin=200)
        protRBFE.inputSetOfMols.set(protExtract)
        protRBFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protRBFE.inputLigand.set('SmallMolecule (g1_1uaz_RET_255-1_1 molecule)')
        protRBFE.ligandB.set('SmallMolecule (g1_1uaz_RET_255-1_1 molecule)')
        protRBFE.setObjLabel('gromacs - pmx RBFE (self-transform sanity check)')

        cls.launchProtocol(protRBFE)
        return protRBFE

    def test(self):
        protImportJNK1 = self._runImportJNK1PDB()
        self._waitOutput(protImportJNK1, 'outputPdb')

        protExtract = self._runExtractLigand(protImportJNK1, jnk1ChainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protRBFE = self._runPmxRBFE(protExtract)
        self._waitOutput(protRBFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protRBFE, 'outputSystem', None))

    def test2(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protRBFE = self._runPmxRBFESelfTransform(protExtract)
        self._waitOutput(protRBFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protRBFE, 'outputSystem', None))


class TestGromacsPmxABFE(TestGromacsPrepareSystem):
    """Smoke test for the pmx absolute binding FEP protocol.

    Tiny/fast window counts and times for testing, not for a physically
    meaningful free energy estimate: real production ABFE needs many more,
    denser windows (especially for the van der Waals decoupling stage) to
    converge."""

    @classmethod
    def _runPmxABFE(cls, protExtract):
        protABFE = cls.newProtocol(
            GromacsPmxABFE, nStepsMin=200, equilTime=2.0, prodTime=4.0,
            nRestrWindows=2, nCoulWindows=2, nVdwWindows=2)
        protABFE.inputSetOfMols.set(protExtract)
        protABFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protABFE.inputLigand.set('SmallMolecule (g1_1uaz_RET_255-1_1 molecule)')
        protABFE.setObjLabel('gromacs - pmx ABFE')

        cls.launchProtocol(protABFE)
        return protABFE

    def test(self):
        protExtract = self._runExtractLigand(self.protImportPDB, chainStr)
        self._waitOutput(protExtract, 'outputSmallMolecules')

        protABFE = self._runPmxABFE(protExtract)
        self._waitOutput(protABFE, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protABFE, 'outputSystem', None))
    test2 = None

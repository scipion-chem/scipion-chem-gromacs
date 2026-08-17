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
from pwchem.protocols.VirtualDrugScreening.protocol_import_smallMolecules import ProtChemImportSmallMolecules
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_preparation import ProtChemOBabelPrepareLigands
from pwchem.protocols.VirtualDrugScreening.protocol_define_manual_structROIs import ProtDefineStructROIs

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

    `test` uses two genuinely different ligands, both docked for real (Vina,
    mirroring pwchem's TestScoreDocking/TestSCORCH2 pipeline: import ->
    OBabel prep -> Vina) into the *same* real binding site of the apo 1uaz
    receptor used to build the GROMACS system - the site is pinned to the
    real centroid of 1uaz chain A's native RET (retinal) ligand (computed
    once from the downloaded 1uaz.cif, see _runDefLigandPocket), not a blind
    whole-protein search. This matters, and was found the hard way: an
    earlier version of this test docked each ligand blindly over the whole
    protein independently, and the two ended up in different sites; pmx's
    atomMapping/ligandHybrid then merged a "unique" atom from ligand B
    straight from its own pose into ligand A's frame with no realignment,
    landing tens of nm away and making mdrun's domain decomposition fail
    outright ("no domain decomposition ... compatible with ... a minimum
    cell size of 47 nm"). RBFE needs co-located poses, not just docked poses.
    Ligands come from the shared "smallMolecules" dataset (4 ZINC compounds).

    Caveat, checked for real (not assumed) with RDKit Morgan-fingerprint
    Tanimoto similarity on the actual mol2 structures: none of the 4 dataset
    ligands are close analogs of one another - all 6 pairwise Tanimoto scores
    fall in 0.057-0.262, below the ~0.4 threshold usually taken as "similar
    enough" for a physically meaningful RBFE alchemical transformation (poor
    forward/reverse work overlap is expected for a low-similarity pair). Of
    the 4, ZINC00000480/ZINC00001019 is the least-dissimilar pair available
    (Tanimoto=0.262), so `test` uses that pair: it exercises the pipeline on a
    real structural change (unlike a self->self transform) without claiming a
    scientifically meaningful RBFE result. See claude/decisions/pmx_RBFE.md.

    `test2` keeps the original self->self (RET->RET) transformation, which
    should give dG(A->B)~=0 within error - a sanity check on the alchemical
    machinery itself, independent of ligand-pair similarity/pose alignment.

    Settings (window counts/times/nStepsMin) are kept tiny/fast for testing
    in both tests, not for a physically meaningful free energy estimate."""

    @classmethod
    def _runDefLigandPocket(cls):
        # Residues 218-226 of chain A: a window around Lys222, the residue covalently
        # bonded (Schiff base) to 1uaz's native RET (retinal) - confirmed for real from
        # the downloaded 1uaz.cif's covale annotations and coordinates (Lys222's real
        # centroid, ~80.7/29.6/5.5, sits a few A from RET's own real centroid,
        # ~90.7/23.6/6.4 - consistent with a direct covalent bond). RET itself is
        # already stripped from protPrepareReceptor's structure (HETATM=True), so this
        # residue window - not the ligand - is what pins both ligands' docking to the
        # same real binding site instead of each drifting off to its own blind pose.
        #
        # surfaceCoords=False: RET sits deeply buried inside 1uaz's transmembrane helix
        # bundle (bacteriorhodopsin), nowhere near the molecular surface - the default
        # surfaceCoords=True tries to snap every input coordinate to the nearest surface
        # point within maxDepth (3 A) and silently yields zero ROIs when none exists
        # nearby, confirmed for real with a single-point "Coordinate:" ROI at RET's own
        # centroid (outputStructROIs came back empty). A multi-atom residue window
        # avoids the surface-mapping step entirely and gives the ROI real spatial
        # extent (needed for a nonzero Vina docking box - a single point has none).
        cls.protPocket = cls.newProtocol(
            ProtDefineStructROIs, inROIs='1) Residues: {"chain": "A", "index": "218-226"}',
            surfaceCoords=False)
        cls.protPocket.inputAtomStruct.set(cls.protPrepareReceptor)
        cls.protPocket.inputAtomStruct.setExtended('outputStructure')
        cls.launchProtocol(cls.protPocket)
        return cls.protPocket

    @classmethod
    def _runImportZincMols(cls):
        dsLig = DataSet.getDataSet('smallMolecules')
        cls.protImportZinc = cls.newProtocol(
            ProtChemImportSmallMolecules, filesPath=dsLig.getFile('mol2'))
        cls.launchProtocol(cls.protImportZinc)
        return cls.protImportZinc

    @classmethod
    def _runZincOBabel(cls, protImport):
        cls.protZincOBabel = cls.newProtocol(
            ProtChemOBabelPrepareLigands, inputType=0, method_charges=0,
            inputSmallMolecules=protImport.outputSmallMolecules, doConformers=False)
        cls.launchProtocol(cls.protZincOBabel)
        return cls.protZincOBabel

    @classmethod
    def _runDockZincMols(cls, protOBabel, protPocket):
        from autodock.protocols import ProtChemVinaDocking
        protVina = cls.newProtocol(
            ProtChemVinaDocking, fromReceptor=1, pocketRadiusN=2, nRuns=1, numberOfThreads=4)
        protVina.inputStructROIs.set(protPocket)
        protVina.inputStructROIs.setExtended('outputStructROIs')
        protVina.inputSmallMolecules.set(protOBabel)
        protVina.inputSmallMolecules.setExtended('outputSmallMolecules')
        cls.launchProtocol(protVina)
        return protVina

    @staticmethod
    def _getMolByZincId(molSet, zincId):
        for mol in molSet:
            if zincId in mol.getMolName():
                return mol.clone()
        raise ValueError(f'{zincId} not found in {molSet}')

    @classmethod
    def _runPmxRBFE(cls, protDock):
        molA = cls._getMolByZincId(protDock.outputSmallMolecules, 'ZINC00000480')
        molB = cls._getMolByZincId(protDock.outputSmallMolecules, 'ZINC00001019')

        protRBFE = cls.newProtocol(
            GromacsPmxRBFE, nStructs=2, swTime=2.0, nvtTime=2.0, nptTime=4.0, nStepsMin=200)
        protRBFE.inputSetOfMols.set(protDock)
        protRBFE.inputSetOfMols.setExtended('outputSmallMolecules')
        protRBFE.inputLigand.set(str(molA))
        protRBFE.ligandB.set(str(molB))
        protRBFE.setObjLabel('gromacs - pmx RBFE (different ligands)')

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
        self._runPrepareReceptor()
        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
        protPocket = self._runDefLigandPocket()
        self._waitOutput(protPocket, 'outputStructROIs', sleepTime=5)

        protImport = self._runImportZincMols()
        self._waitOutput(protImport, 'outputSmallMolecules')
        protOBabel = self._runZincOBabel(protImport)
        self._waitOutput(protOBabel, 'outputSmallMolecules')
        protDock = self._runDockZincMols(protOBabel, protPocket)
        self._waitOutput(protDock, 'outputSmallMolecules', sleepTime=10)

        protRBFE = self._runPmxRBFE(protDock)
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

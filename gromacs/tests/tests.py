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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

from gromacs.protocols import GromacsSystemPrep, GromacsModifySystem, GromacsMDSimulation


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
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_prepared.pdb'))
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPrepareSystem(cls):
        protPrepare = cls.newProtocol(
            GromacsSystemPrep,
            inputStructure=cls.protImportPDB.outputPdb,
            boxType=1, sizeType=1, padDist=5.0,
            mainForceField=0, waterForceField=2,
            placeIons=1, cationType=7, anionType=1)

        cls.launchProtocol(protPrepare)
        return protPrepare

    def test(self):
        protPrepare = self._runPrepareSystem()
        self._waitOutput(protPrepare, 'outputSystem', sleepTime=10)
        self.assertIsNotNone(getattr(protPrepare, 'outputSystem', None))


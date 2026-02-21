# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     James Krieger (jmkrieger@cnb.csic.es)
# *
# * Centro Nacional de Biotecnologia, CSIC
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
This module will run the EMMIVox workflow to refine single structures or ensembles
against cryo-EM maps using gromacs and plumed.
"""
from os.path import basename, splitext, join

from pwem.convert.atom_struct import cifToPdb
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol

from pyworkflow.protocol import params

from gromacs import Plugin
from gromacs.objects import GromacsSystem

SUMMARY_NO_OUTPUT = 'Output structure not ready yet'

class EMMIVoxProtocol(EMProtocol):
    """
    This module will run the EMMIVox workflow to refine single structures or ensembles
    against cryo-EM maps using gromacs and plumed.
    """
    _label = 'EMMIVox'
    _possibleOutputs = {'outputStructure': AtomStruct, 'outputSystem': GromacsSystem}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label='EMMIVox options')
        
        form.addParam('inputStructure', params.PointerParam, label="Reference structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='This should be an AtomStruct corresponding to your starting structure that fits the map '
                            'and has the right residue numbers.')
        
        form.addParam('inputSystem', params.PointerParam, label="Prepared MD system",
                      important=True,
                      pointerClass='GromacsSystem',
                      help='This should be based on the reference structure, but it may have moved in space and/or been renumbered.')

        form.addParam('inputVolumes', params.PointerParam, label="Target volumes",
                      important=True, allowsNull=True,
                      pointerClass='Volume',
                      help='This should be a main map and two half maps that the reference structure is already fitted to')

        form.addParam('cutoff', params.FloatParam, label='Correlated voxel cutoff',
                      important=True, default=0.9,
                      help='This cutoff will be used to exclude correlated voxels (above this threshold). '
                      'If you want to keep all the voxels of the input map, set this to 1.0')

        form.addParam('zoneMap', params.BooleanParam, label='Zone volume?',
                      important=True, default=True,
                      help='This will zone the volume close to the reference model (optional, but speeds up things a lot)')

        form.addParam('zoneDistance', params.FloatParam, label='Volume zoning distance',
                      important=True, default=3.5,
                      help='This distance cutoff will be used to exclude correlated voxels (above this threshold). '
                      'If you want to keep all the voxels of the input map, set this to 1.0')

        form.addParam('selection1', params.StringParam, default="protein",
                      label="Selection string for atoms in map",
                      help='MDAnalysis selection string for the atoms from the MD system that will be used to simulate the volume for comparison. '
                            'This will also be used to align the reference model and volumes to the MD system')

        form.addParam('selection2', params.StringParam, default="protein",
                      label="Selection string for reference structure",
                      help='MDAnalysis selection string for the atoms in the reference structure to align it and the volumes to the MD system')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        inputStructFn = self.inputStructure.get().getFileName()
        if inputStructFn.endswith('.cif'):
            self.inputPdbFn = inputStructFn[:-4]+'.pdb'    
            cifToPdb(inputStructFn, self.inputPdbFn)
        else:
            self.inputPdbFn = inputStructFn

        inputSystem = self.inputSystem.get()
        if inputSystem.hasTrajectory() and self.prevTrj.get():
            self.inputGroFn = inputSystem.getOriStructFile()
        else:
            self.inputGroFn = inputSystem.getSystemFile()
        self.inputTopFn = inputSystem.getTopologyFile()

        # step 0, system setup, stage 1: renumber gro and pdb files using topology
        self.outputPdbFn = self._getExtraPath(splitext(basename(self.inputPdbFn))[0] + '.pdb')
        self.outputGroFn = self._getExtraPath(splitext(basename(self.inputGroFn))[0] + '.gro')

        tutorialLocation = Plugin._getEmmiVoxTutorialLocation()
        STEP1_0_LOCATION = join(tutorialLocation, '1-refinement', '0-Building')

        args = '{0} {1} {2} {3} {4}'.format(self.inputGroFn, self.inputPdbFn, STEP1_0_LOCATION,
                                            self._getExtraPath(), self.inputTopFn)
        self.runJob(Plugin.getEmmiVoxProgram('renumber.sh',
                                             location='tutorials/1-refinement/0-Building'),
                    args)

        # step 0, system setup, stage 2: create index groups for simulating the map and writing to xtc
        self.runJob(Plugin.getEmmiVoxProgram('make_ndx.py', python=True),
                    '{0} {1} System-MAP --ndx index.ndx'.format(self.inputGroFn,
                                                                self.selection1.get()))
        self.runJob(Plugin.getEmmiVoxProgram('make_XTC_ndx.py', python=True),
                    'System-MAP-H --ndx index.ndx')

        # step 1, stage 1: map preparation
        volumes = self.inputVolumes.get()
        mainMap = volumes.getFileName()
        halfMaps = volumes.getHalfMaps().split(',')
        refPdbFn = self.inputStructure.get().getFileName()
        selstr1 = self.selection1.get()
        selstr2 = self.selection2.get()

        args = '{0} {1} emd_plumed.dat --halfmaps {2}'.format(mainMap, self.cutoff.get(),
                                                              ' '.join(halfMaps))
        if self.zoneMap.get():
            args += ' --zone {0} --zone_PDB {1} --zone_sel {2}'.format(self.zoneDistance.get(),
                                                                      refPdbFn,
                                                                      selstr1)
        self.runJob(Plugin.getEmmiVoxProgram('cryo-EM_preprocess.py', python=True), args)

        # step 1, stage 2: map and model alignment
        args = '{0} {1} emd_plumed_aligned.dat emd_plumed.dat --ref_sel {2} --mobile_sel {3}'.format(self.outputPdbFn,
                                                                                                     refPdbFn, selstr1,
                                                                                                     selstr2)
        self.runJob(Plugin.getEmmiVoxProgram('align-VOXELS.py', python=True), args)

    def createOutputStep(self):
        outputPdb = AtomStruct()
        outputPdb.setFileName(self.outputPdbFn)
        self._defineOutputs(outputStructure=outputPdb)

        outSystem = GromacsSystem(filename=self.outputGroFn,
                                  oriStructFile=self.outputGroFn,
                                  topoFile=self.inputTopFn)
        self._defineOutputs(outputSystem=outSystem)

    def _summary(self):
        if not hasattr(self, 'outputStructure'):
            summ = [SUMMARY_NO_OUTPUT]
        else:
            summ = ['The new structure has been created']
        return summ

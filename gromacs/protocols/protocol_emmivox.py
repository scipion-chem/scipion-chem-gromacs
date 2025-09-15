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

from pyworkflow.protocol.params import PointerParam

from gromacs.constants import *
from gromacs import Plugin
from gromacs.objects import *

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
        
        form.addParam('inputStructure', PointerParam, label="Reference structure",
                      important=True,
                      pointerClass='AtomStruct',
                      help='This should be an AtomStruct corresponding to your starting structure that fits the map '
                            'and has the right residue numbers.')
        
        form.addParam('inputSystem', PointerParam, label="Prepared MD system",
                      important=True,
                      pointerClass='GromacsSystem',
                      help='This should be based on the reference structure, but it may have moved in space and/or been renumbered.')

        form.addParam('inputVolumes', PointerParam, label="Target volumes",
                      important=True, allowsNull=True,
                      pointerClass='Volume',
                      help='This should be a main map and two half maps that the reference structure is already fitted to')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('computeStep')
        self._insertFunctionStep('createOutputStep')

    def computeStep(self):
        inputStructFn = self.inputStructure.get().getFileName()
        if inputStructFn.endswith('.cif'):
            self.inputPdbFn = inputStructFn[:-4]+'.pdb'    
            cifToPdb(inputStructFn, self.inputPdbFn)
        else:
            self.inputPdbFn = inputStructFn

        inputSystem = self.inputSystem.get()
        if self.gromacsSystem.get().hasTrajectory() and self.prevTrj.get():
            self.inputGroFn = self.gromacsSystem.get().getOriStructFile()
        else:
            self.inputGroFn = self.gromacsSystem.get().getSystemFile()
        self.inputTopFn = inputSystem.getTopologyFile()

        self.outputPdbFn = self._getPath(splitext(basename(self.inputPdbFn))[0] + '_fixed.pdb')
        self.outputGroFn = self._getPath(splitext(basename(self.inputGroFn))[0] + '_fixed.gro')

        tutorialLocation = Plugin._getEmmiVoxTutorialLocation()
        step1_0_location = join(tutorialLocation, '1-refinement', '0-Building')

        args = '{0} {1} {2} {3} {4}'.format(self.inputGroFn, self.inputPdbFn, step1_0_location,
                                        self._getExtraPath(), self.inputTopFn)
        self.runJob(Plugin.getEmmiVoxProgram('renumber.sh',
                                             location='tutorials/1-refinement/0-Building'),
                    args)

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

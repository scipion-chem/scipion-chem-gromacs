# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *          	 James M. Krieger (jamesmkrieger@gmail.com)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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

DISTANCE = 0
ANGLE = 1
DIHEDRAL = 2

selstrHelp = '''The distance, angle, or dihedral will be calculated using the centers of the selections.
Selections are defined using comma-separated lists of ProDy selection strings
(see http://http://www.bahargroup.org/prody/tutorials/prody_tutorial/selection.html)
but are reinterpreted using MDTraj and Biopython to get atom IDs.
Each of the comma-separated selection strings can only have one chain name.'''

defaultSelstr1 = "chain A and name CA and resid 117 to 243,chain A and name CA and resid 354 to 380"
defaultSelstr2 = "chain A and name CA and resid 4 to 116,chain A and name CA and resid 244 to 353"
defaultSelstr3 = "chain B and name CA and resid 4 to 116,chain B and name CA and resid 244 to 353"
defaultSelstr4 = "chain B and name CA and resid 117 to 243,chain B and name CA and resid 354 to 380"
measureTypeCheck = "measureType>%d"

GROUP_ATOMS_STR = '{0}: CENTER ATOMS={1}\n'
PRINT_STR = 'PRINT STRIDE=1 ARG={0} FILE=COLVAR FMT=%6.3f\n'

"""
This module will perform analyses using plumed programs, e.g. driver and sum_hills
"""
from Bio.PDB.PDBParser import PDBParser
import mdtraj
import numpy as np

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwem.convert.atom_struct import cifToPdb
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol

from gromacs.objects import *
from gromacs.constants import *
from gromacs import Plugin as gromacsPlugin

plumed_ver = gromacsPlugin._getInstalledVersion(PLUMED_DIC)

if plumed_ver.startswith('2'):
	PLUMED_DOC_LINK = 'https://www.plumed.org/doc-v{0}/user-doc/html/_g_h_o_s_t.html'.format(plumed_ver[:-1])
else:
	PLUMED_DOC_LINK = 'https://www.plumed.org/doc-master/user-doc/html/GHOST/'

class PlumedRunAnalysis(EMProtocol):
	"""
	This protocol will analyse structures and simulations with plumed.
	"""
	_label = 'Plumed Analysis'

	_analysisTypes = ['driver']#, 'sum_hills']

	# -------------------------- DEFINE constants ----------------------------
	def __init__(self, **kwargs):
		EMProtocol.__init__(self, **kwargs)


	# -------------------------- DEFINE param functions ----------------------
	def _defineParams(self, form):

		""" Define the input parameters that will be used.
		"""
		form.addSection(label=Message.LABEL_INPUT)

		form.addParam('inputSystem', params.PointerParam, label="Input MD System: ",
					  pointerClass='MDSystem,AtomStruct',
					  help='MD system to be analysed')

		group = form.addGroup('Analysis')
		group.addParam('analysisType', params.EnumParam,
					   label='Analysis type: ',
					   choices=self._analysisTypes, default=0,
					   help='Type of analysis to perform in a step')

		form.addParam('measureType', params.EnumParam,
					  choices=['distance', 'angle', 'dihedral'], default=DISTANCE,
					  label='Measure type',
					  help='Select the type of measure for defining a collective variable.')

		group = form.addGroup('Selection 1')
		group.addParam('selection1', params.StringParam, default=defaultSelstr1,
					  label="Selection string 1",
					  help=selstrHelp)
		group.addParam('label1', params.StringParam, default='com1',
					  label="Label for selection 1 center",
					  help=labelHelp)

		group = form.addGroup('Selection 2')
		group.addParam('selection2', params.StringParam, default=defaultSelstr2,
					  label="Selection string 2",
					  help=selstrHelp)
		group.addParam('label2', params.StringParam, default='com2',
					  label="Label for selection 2 center",
					  help=labelHelp)
		
		group = form.addGroup('Selection 3', condition=measureTypeCheck % DISTANCE)
		group.addParam('selection3', params.StringParam, default=defaultSelstr3,
					  label="Selection string 3", condition=measureTypeCheck % DISTANCE,
					  help=selstrHelp)
		group.addParam('label3', params.StringParam, default='com3',
					  label="Label for selection 3 center", condition=measureTypeCheck % DISTANCE,
					  help=labelHelp)

		group = form.addGroup('Selection 4', condition=measureTypeCheck % ANGLE)
		group.addParam('selection4', params.StringParam, default=defaultSelstr4,
					  label="Selection string 4", condition=measureTypeCheck % ANGLE,
					  help=selstrHelp)
		group.addParam('label4', params.StringParam, default='com4',
					  label="Label for selection 4 center", condition=measureTypeCheck % ANGLE,
					  help=labelHelp)

	# --------------------------- STEPS functions ------------------------------
	def _insertAllSteps(self):
		# Insert processing steps
		self._insertFunctionStep('analyseStep')
		self._insertFunctionStep('createOutputStep')

	def analyseStep(self):
		self.plumedFile = self.generatePlumedFile()
		progType = self._analysisTypes[self.analysisType.get()]
		self.callPlumed(self.plumedFile, progType)

	def createOutputStep(self):
		inputSystem = self.inputSystem.get()
		if isinstance(inputSystem, MDSystem):
			outSystem = GromacsSystem(filename=inputSystem.getFileName(),
									topoFile=inputSystem.getTopologyFile(),
									colvarFile=self.getColvarFile(),
									trjFile=inputSystem.getTrajectoryFile())
			self._defineOutputs(outputSystem=outSystem)
		else:
			outStructure = AtomStruct(filename=self.getPdbFile())
			self._defineOutputs(outputSystem=outStructure)

	def generatePlumedFile(self):
		plumedFile = self._getPath('plumed.dat')
		plumedStr = ''

		idLists = self.getChainNames()
		mdtrajTop = mdtraj.load(self.getPdbFile()).top

		measureType = self.measureType.get()
		atomIds1 = self.selectAtomIds(self.selection1.get(), mdtrajTop, idLists)
		plumedStr += GROUP_ATOMS_STR.format(self.label1.get(), atomIds1)
		
		atomIds2 = self.selectAtomIds(self.selection2.get(), mdtrajTop, idLists)
		plumedStr += GROUP_ATOMS_STR.format(self.label2.get(), atomIds2)

		if measureType == DISTANCE:
			plumedStr += 'd0: DISTANCE ATOMS={0},{1}\n'.format(self.label1.get(), self.label2.get())
			plumedStr += PRINT_STR.format('d0')
		else:
			atomIds3 = self.selectAtomIds(self.selection3.get(), mdtrajTop, idLists)
			plumedStr += GROUP_ATOMS_STR.format(self.label3.get(), atomIds3)

			if measureType == ANGLE:
				plumedStr += 'a0: ANGLE ATOMS={0},{1},{2}\n'.format(self.label1.get(), self.label2.get(),
													self.label3.get())
				plumedStr += PRINT_STR.format('a0')
			else:
				atomIds4 = self.selectAtomIds(self.selection4.get(), mdtrajTop, idLists)
				plumedStr += GROUP_ATOMS_STR.format(self.label4.get(), atomIds4)

				if measureType == DIHEDRAL:
					plumedStr += 't0: TORSION ATOMS={0},{1},{2},{3}\n'.format(self.label1.get(), self.label2.get(),
														self.label3.get(), self.label4.get())
					plumedStr += PRINT_STR.format('t0')

		with open(plumedFile, 'w') as f:
			f.write(plumedStr)
		return plumedFile

	def generatePdbFile(self):
		inFileName = os.path.abspath(self.inputSystem.get().getFileName())
		pdbFileName = self._getPath(os.path.basename(inFileName))[:-4] + '.pdb'
		if inFileName.endswith('.cif'):
			cifToPdb(inFileName, pdbFileName)
		else:
			args = ' trjconv -s {} -f {} -o {}'.format(inFileName, inFileName,
											  pdbFileName)
			gromacsPlugin.runGromacsPrintf(printfValues=['System'], args=args,
								  cwd=self._getPath())
		return pdbFileName

	def callPlumed(self, plumedFile, progType):
		"""Prepare input files and call plumed"""
		if isinstance(self.inputSystem.get(), MDSystem):
			trjFile = self.inputSystem.get().getTrajectoryFile()
			if not trjFile:
				trjFile = self.inputSystem.get().getSystemFile()
		else:
			trjFile = self.inputSystem.get().getFileName()

		command = '%s --plumed %s' % (progType, os.path.abspath(plumedFile))

		trjType = os.path.splitext(trjFile)[1][1:]

		if trjType == 'cif':
			trjFile = self.getPdbFile()
			trjType = 'pdb'

		if trjType in ['xyz', 'gro', 'dlp4', 'xtc', 'trr']:
			trjTypeFlag = 'i' + trjType
		else:
			trjTypeFlag = 'mf_' + trjType

		command += ' --%s %s' % (trjTypeFlag, os.path.abspath(trjFile))
		trjFilePath = os.path.split(trjFile)[0]
		gromacsPlugin.runPlumed(self, 'plumed', command, cwd=trjFilePath, numberOfMpi=0)

		shutil.move(os.path.join(trjFilePath, 'COLVAR'), self.getColvarFile())

	def getColvarFile(self):
		return self._getPath('COLVAR')

	def getPdbFile(self):
		inFileName = self.inputSystem.get().getFileName()
		pdbFileName = self._getPath(os.path.basename(inFileName))[:-4] + '.pdb'

		if os.path.exists(pdbFileName):
			return pdbFileName
		
		if inFileName.endswith('.pdb'):
			shutil.copy(inFileName, pdbFileName)
			return pdbFileName
		
		if os.path.splitext(inFileName)[1] in ['.cif', '.gro']:
			pdbFileName = self.generatePdbFile()

		return pdbFileName

	def getChainNames(self):
		inputSystem = self.inputSystem.get()
		if hasattr(inputSystem, 'getChainNames') and inputSystem.getChainNames():
			return inputSystem.getChainNames()

		parser = PDBParser()
		structure = parser.get_structure('', self.getPdbFile())
		chainNames = [ch.get_id() for ch in list(structure.get_models())[0].get_chains()]

		resids = [''.join([str(x) for x in r.get_id()]).strip() for r in structure.get_residues()]
		atomids = [a.get_serial_number() for a in structure.get_atoms()]

		return chainNames, resids, atomids

	def selectAtomIds(self, selstrs, mdtrajTop, idLists):
		chainNames, resids, atomids = idLists

		atomIds1 = []
		for selstr in [selstr.strip() for selstr in selstrs.split(',')]:
			chainStrStart = selstr.find('chain')
			chainStrEnd = selstr.find('chain') + len('chain ')
			chainNameStr = selstr[chainStrEnd: chainStrEnd + selstr[chainStrEnd:].find(' ')]
			chainId = chainNames.index(chainNameStr)
			chainIdSel = 'chainid {0}'.format(chainId)
			
			selstr = selstr.replace(selstr[chainStrStart:chainStrEnd+len(chainNameStr)], chainIdSel)

			selstr = selstr.replace('resnum','resid').replace(' to ','-')
			
			residsStrEnd = selstr.find('resid') + len('resid ')
			if selstr[residsStrEnd:].find(' ') != -1:
				residNameStr = selstr[residsStrEnd: residsStrEnd + selstr[residsStrEnd:].find(' ')]
			else:
				residNameStr = selstr[residsStrEnd:]

			selstr = selstr[:residsStrEnd] + ' to '.join([str(np.nonzero(np.array(resids) == residName)[0][chainId])
												 		  for residName in residNameStr.split('-')])

			atomIds1.extend(list(np.array(atomids)[mdtrajTop.select(selstr)]))

		return ','.join([str(atomId) for atomId in atomIds1])

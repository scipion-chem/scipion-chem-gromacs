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
This module will prepare the system for the simulation
"""
import os, subprocess, shutil, uuid, json

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.convert import AtomicStructHandler

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import runOpenBabel, removeStartWithLines, cleanPDB
from pwchem.protocols import ProtocolLigandParametrization
from pwchem.constants import OPENBABEL_DIC, RDKIT_DIC

from gromacs import Plugin as gromacsPlugin
import gromacs.objects as grobj
from gromacs.constants import *
from Bio import PDB

GROMACS_AMBER03 = 0
GROMACS_AMBER94 = 1
GROMACS_AMBER96 = 2
GROMACS_AMBER99 = 3
GROMACS_AMBER99SB = 4
GROMACS_AMBERSB_ILDN = 5
GROMACS_AMBERGS = 6
GROMACS_CHARMM27 = 7
GROMACS_CHARMM36 = 8
GROMACS_GROMOS43A1 = 9
GROMACS_GROMOS43A2 = 10
GROMACS_GROMOS45A3 = 11
GROMACS_GROMOS53A5 = 12
GROMACS_GROMOS53A6 = 13
GROMACS_GROMOS54A7 = 14
GROMACS_OPLSAA = 15
GROMACS_AMBER14SB = 16
GROMACS_AMBER19SB = 17

GROMACS_MAINFF_NAME = dict()
GROMACS_MAINFF_NAME[GROMACS_AMBER03] = 'amber03'
GROMACS_MAINFF_NAME[GROMACS_AMBER94] = 'amber94'
GROMACS_MAINFF_NAME[GROMACS_AMBER96] = 'amber96'
GROMACS_MAINFF_NAME[GROMACS_AMBER99] = 'amber99'
GROMACS_MAINFF_NAME[GROMACS_AMBER99SB] = 'amber99sb'
GROMACS_MAINFF_NAME[GROMACS_AMBERSB_ILDN] = 'amber99sb-ildn'
GROMACS_MAINFF_NAME[GROMACS_AMBERGS] = 'amberGS'
GROMACS_MAINFF_NAME[GROMACS_CHARMM27] = 'charmm27'
GROMACS_MAINFF_NAME[GROMACS_CHARMM36] = 'charmm36-feb2021'
GROMACS_MAINFF_NAME[GROMACS_GROMOS43A1] = 'gromos43a1'
GROMACS_MAINFF_NAME[GROMACS_GROMOS43A2] = 'gromos43a2'
GROMACS_MAINFF_NAME[GROMACS_GROMOS45A3] = 'gromos45a3'
GROMACS_MAINFF_NAME[GROMACS_GROMOS53A5] = 'gromos53a5'
GROMACS_MAINFF_NAME[GROMACS_GROMOS53A6] = 'gromos53a6'
GROMACS_MAINFF_NAME[GROMACS_GROMOS54A7] = 'gromos54a7'
GROMACS_MAINFF_NAME[GROMACS_OPLSAA] = 'oplsaa'
GROMACS_MAINFF_NAME[GROMACS_AMBER14SB] = 'amber14sb'
GROMACS_MAINFF_NAME[GROMACS_AMBER19SB] = 'amber19sb'

GROMACS_LIST = [GROMACS_MAINFF_NAME[GROMACS_AMBER03], GROMACS_MAINFF_NAME[GROMACS_AMBER94],
                GROMACS_MAINFF_NAME[GROMACS_AMBER96], GROMACS_MAINFF_NAME[GROMACS_AMBER99],
                GROMACS_MAINFF_NAME[GROMACS_AMBER99SB], GROMACS_MAINFF_NAME[GROMACS_AMBERSB_ILDN],
                GROMACS_MAINFF_NAME[GROMACS_AMBERGS],
                GROMACS_MAINFF_NAME[GROMACS_CHARMM27], GROMACS_MAINFF_NAME[GROMACS_CHARMM36],
                GROMACS_MAINFF_NAME[GROMACS_GROMOS43A1], GROMACS_MAINFF_NAME[GROMACS_GROMOS43A2],
                GROMACS_MAINFF_NAME[GROMACS_GROMOS45A3], GROMACS_MAINFF_NAME[GROMACS_GROMOS53A5],
                GROMACS_MAINFF_NAME[GROMACS_GROMOS53A6], GROMACS_MAINFF_NAME[GROMACS_GROMOS54A7],
                GROMACS_MAINFF_NAME[GROMACS_OPLSAA],
                GROMACS_MAINFF_NAME[GROMACS_AMBER14SB], GROMACS_MAINFF_NAME[GROMACS_AMBER19SB]]

NUCLEIC_RESNAMES = {'A', 'C', 'G', 'U', 'T', 'DA', 'DC', 'DG', 'DT', 'DU'}

GROMACS_SPC = 0
GROMACS_SPCE = 1
GROMACS_TIP3P = 2
GROMACS_TIP4P = 3
GROMACS_TIP5P = 4
GROMACS_OPC = 5
GROMACS_OPC3 = 6

GROMACS_WATERFF_NAME = dict()
GROMACS_WATERFF_NAME[GROMACS_SPC] = 'spc'
GROMACS_WATERFF_NAME[GROMACS_SPCE] = 'spce'
GROMACS_WATERFF_NAME[GROMACS_TIP3P] = 'tip3p'
GROMACS_WATERFF_NAME[GROMACS_TIP4P] = 'tip4p'
GROMACS_WATERFF_NAME[GROMACS_TIP5P] = 'tip5p'
GROMACS_WATERFF_NAME[GROMACS_OPC] = 'opc'
GROMACS_WATERFF_NAME[GROMACS_OPC3] = 'opc3'

GROMACS_WATERS_LIST = [GROMACS_WATERFF_NAME[GROMACS_SPC], GROMACS_WATERFF_NAME[GROMACS_SPCE],
GROMACS_WATERFF_NAME[GROMACS_TIP3P], GROMACS_WATERFF_NAME[GROMACS_TIP4P],
GROMACS_WATERFF_NAME[GROMACS_TIP5P],
GROMACS_WATERFF_NAME[GROMACS_OPC], GROMACS_WATERFF_NAME[GROMACS_OPC3]]

WATER_BOXES = {'spc': 'spc216', 'spce': 'spc216', 'tip3p': 'spc216', 'opc3': 'spc216',
               'tip4p': 'tip4p', 'tip4pew': 'tip4p', 'opc': 'tip4p', 'tip5p': 'tip5p'}

STRUCTURE, LIGAND = 0, 1
TOPOL_TOP = "topol.top"

GAPS_OPTIONS = ['No', 'Gaps termini', 'All termini']
SSBONDS_OPTIONS = ['None', 'Automatic', 'Manual']

def getMol2FormalCharge(mol2File):
  '''Net formal charge of a mol2 file. OpenBabel writes the formal charges in @<TRIPOS>UNITY_ATOM_ATTR, a block
  Antechamber ignores, so without passing it ACPYPE guesses 0 even for charged ligands'''
  charge, inAttrs = 0, False
  with open(mol2File) as f:
    for line in f:
      if line.startswith('@<TRIPOS>'):
        inAttrs = line.strip() == '@<TRIPOS>UNITY_ATOM_ATTR'
      elif inAttrs:
        fields = line.split()
        if len(fields) == 2 and fields[0] == 'charge':
          charge += int(fields[1])
  return charge

def replaceInFile(file, inStr, repStr):
  inStr, repStr = inStr.replace('\n', '\\n'), repStr.replace('\n', '\\n')
  subprocess.check_call(f"sed -i -z 's/{inStr}/{repStr}/g' {file}", shell=True)
  return file


class GromacsSystemPrep(ProtocolLigandParametrization):
    """
    AI Generated:

    Protocol to prepare a GROMACS system ready for Molecular Dynamics.

    It builds the topology, coordinates, simulation box, solvent and ions of a protein
    (optionally in complex with a docked ligand) with the GROMACS tools, producing a
    GromacsSystem that can be directly simulated with "Run MD simulation".

    Overview
    --------
    The protocol performs the following steps:

    1. Optional structure preprocessing
       - PDBFixer: adds missing atoms and residues, replaces non-standard residues
       - Removes waters, hetero residues (e.g. a crystal ligand) and unselected chains

    2. Optional ligand parametrization
       - The selected docked ligand is parametrized with ACPYPE (GAFF/GAFF2)

    3. Topology generation with pdb2gmx
       - Applies the chosen protein and water force fields
       - Optionally adds ACE/NME caps to the gaps and/or termini
       - Optionally defines the disulfide bonds

    4. Box definition, solvation and ions addition
       - editconf, solvate and genion

    5. Output creation
       - GromacsSystem with topology, coordinates, restraints and index file

    Input
    -----
    inputFrom:
        Type of input:
        - AtomStruct: a protein structure (apo system)
        - SetOfSmallMolecules: a set of docked molecules; the receptor is taken from the
          set and one ligand (inputLigand) is prepared in complex with it

    inputStructure / inputSetOfMols, inputLigand:
        Structure, or set of docked molecules and the name of the ligand to prepare.

    Preprocess structure (waters, HETATM, het2keep, rchains, chain_name):
        Same cleaning as pwchem's target preparation. Removing the HETATM is needed when the
        receptor still holds a crystal ligand: pdb2gmx has no topology for it. Ions or cofactors
        to simulate can be kept with het2keep.

    usePDBFixer (addAtoms, addRes, repNonStd, extraClean):
        Repair the input structure before the cleaning. "Heavy" atoms is recommended,
        since pdb2gmx adds the hydrogens of the chosen force field itself. Missing residues
        are only detected from the SEQRES records of the input.

    mainForceField / waterForceField:
        Protein force field (AMBER, CHARMM, GROMOS, OPLS-AA families) and water model.

    addCaps:
        Add ACE/NME capping groups to the internal gaps or to all the termini.

    Ligand parametrization (chargeMethod, netCharge, multip, atomType, qprog):
        ACPYPE options. If netCharge is empty, the formal charge written in the ligand file
        is used.

    Boundary box (boxType, sizeType, distances):
        Cubic or orthorhombic box, defined by its absolute size or by the buffer distance
        between the solute and the box edges.

    Ions (placeIons, cationType, anionType, numbers, addSalt, saltConc):
        Neutralize the system, add a fixed number of ions and/or a salt concentration.

    SS bonds (handleSSBonds, selectSSBonds):
        Define none, all the detected or a manual selection of disulfide bonds.

    Workflow
    --------
    1. preprocessStructureStep (optional):
       - Runs PDBFixer and/or removes waters, hetero residues and chains from the input receptor;
         the resulting file is used in every later step

    2. parametrizeLigandStep (only with docked molecules):
       - Adds hydrogens to the ligand and runs ACPYPE

    3. PDB2GMXStep:
       - Caps the protein if requested and renumbers the residues if several chains
         share residue numbers
       - Runs pdb2gmx (retrying with -ignh if the input hydrogens are not recognized)
       - Merges the ligand topology and coordinates into the system

    4. editConfStep / solvateStep / addIonsStep:
       - Defines the box, solvates it with the water model and adds the ions

    5. createOutputStep:
       - Creates the GromacsSystem and its index file, including the protein-ligand
         complex group (Protein_<ligand>) named after the ligand residue

    Output
    ------
    outputSystem:
        GromacsSystem containing:
        - solvated coordinates (.gro) and topology (.top)
        - position restraints file and index file (.ndx)
        - chain names and lengths, force field and water model
        - ligand topology and residue name (when a ligand is prepared)

    Key Features
    ------------
    - Apo systems and protein-ligand complexes from docking results
    - Structure cleaning and repair with PDBFixer
    - Wide choice of GROMACS force fields and water models
    - Capping, disulfide bond and ion handling

    Notes
    -----
    - The form checks that the water model, the ions to add and the ions present in the
      input are supported by the chosen force field, and warns when the ion parameters
      were fitted for a different water model
    - Non-standard residues or ligands not passed as docked molecules cannot be
      parametrized by pdb2gmx
    """
    _label = 'System preparation'
    IMPORT_FROM_FILE = 0
    IMPORT_FROM_SCIPION = 1
    IMPORT_MDP_FILE = 0
    IMPORT_MDP_SCIPION = 1

    _cations = [CA, CS, CU, CU2, K, LI, MG, NA, RB, ZN]
    _anions = [BR, CL, F, I]

    # -------------------------- DEFINE constants ----------------------------


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('inputFrom', params.EnumParam, default=STRUCTURE, display=params.EnumParam.DISPLAY_HLIST,
                      label='Input from: ', choices=['AtomStruct', 'SetOfSmallMolecules'],
                      help='Type of input you want to use')
        form.addParam('inputStructure', params.PointerParam, pointerClass='AtomStruct',
                      label='Input structure to be prepared for MD:', allowsNull=False, condition='inputFrom==0',
                      help='Atomic structure to be prepared for MD by solvation, ions addition etc')
        form.addParam('inputSetOfMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label='Input set of molecules:', allowsNull=False, condition='inputFrom==1',
                      help='Input set of docked molecules. One of them will be prepared together with its target')
        form.addParam('inputLigand', params.StringParam, condition='inputFrom==1',
                      label='Ligand to prepare: ',
                      help='Specific ligand to prepare in the system')

        group = form.addGroup('Preprocess structure')
        self._definePreprocessParams(group)

        group = form.addGroup('Force field')
        self._defineFFParams(group)

        form.addParam('addCaps', params.EnumParam, choices=GAPS_OPTIONS, default=0,
                       label='Add ACE and NME caps: ',
                       help='Add acetyl (ACE) and N-methylamide (NME) capping groups to protein N-termini and C-termini respectively. '
                            'These caps neutralize terminal charges and are commonly used in MD simulations. '
                            '\n*None*: No caps added. '
                            '\n*Gaps termini*: Add caps only to missing loops (internal gaps in the structure), '
                            'preserving the real N- and C-termini uncapped. '
                            '\n*All termini*: Add caps to both gaps and real protein termini.')

        self._defineACPYPEparams(form, condition=f'inputFrom=={LIGAND}')
        form.getParam('netCharge').help.set(
          'Net molecular charge (int). Leave empty (recommended) to use the formal charge written in the ligand file. '
          'Set it only to override that charge, e.g. if the input file has no charge information (pdbqt) or to use a '
          'different protonation state.')

        form.addSection('MD prep')
        group = form.addGroup('Boundary box')
        self._defineBoxParams(group)

        group = form.addGroup('Ions')
        self._defineIonsParams(group)

        group = form.addGroup('SS bonds')
        self._defineSSBondsParams(group)

    def _definePreprocessParams(self, group):
        """Cleaning options of pwchem's target preparation (same names and defaults) followed by PDBFixer"""
        group.addParam("waters", params.BooleanParam, label='Remove waters: ', default=False,
                       help='Remove all water molecules from the input structure (the system is solvated later)')
        group.addParam("HETATM", params.BooleanParam, label='Remove ligands HETATM: ', default=False,
                       help='Remove the ligands and HETATM contained in the input structure (e.g. a co-crystallized '
                            'ligand), as pdb2gmx can only parametrize the residues of the force field.\n'
                            'Ions and cofactors to simulate can be kept with the next parameter.')
        group.addParam('het2keep', params.StringParam, condition='HETATM', default='',
                       label='Heteroatom residues to keep: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Heteroatom residues to keep in the structure (comma-separated), e.g. metal ions such '
                            'as ZN or MG. They must be parametrized by the chosen force field.')
        group.addParam("rchains", params.BooleanParam, label='Select specific chains: ', default=False,
                       help='Keep only the chains selected')
        group.addParam("chain_name", params.StringParam, label="Keep chains: ", condition="rchains==True",
                       help="Select the chain(s) you want to keep in the structure")
        self._definePDBFixerParams(group)

    def _definePDBFixerParams(self, group):
        group.addParam("usePDBFixer", params.BooleanParam, label='Use PDBFixer: ', default=False,
                       help='Repair the input structure with PDBFixer before the system preparation: add missing '
                            'atoms and residues and replace non-standard residues.')
        group.addParam('addAtoms', params.EnumParam, default=1, condition="usePDBFixer",
                       label="Add missing atoms: ", choices=['All', 'Heavy', 'Hydrogen', 'None'],
                       help='Missing atoms to add with PDBFixer.\nHeavy is recommended: pdb2gmx adds the hydrogens '
                            'of the chosen force field itself, and hydrogens named by PDBFixer may not match it.')
        group.addParam('addRes', params.BooleanParam, default=False,
                       label="Add missing residues: ", condition="usePDBFixer",
                       help='Use PDBFixer to add missing residues. They are only detected from the SEQRES records '
                            'of the input PDB, so structures without them (e.g. predicted models) are not modified.')
        group.addParam('repNonStd', params.BooleanParam, default=False,
                       label="Replace non-standard residues: ", condition="usePDBFixer",
                       help='Use PDBFixer to replace nonstandard residues with standard equivalents, '
                            'which pdb2gmx would not be able to parametrize')
        group.addParam('extraClean', params.BooleanParam, default=False, condition="usePDBFixer",
                       label='Remove HETATM after PDBFixer: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Remove the HETATM records after PDBFixer has been used (including ions and waters). '
                            'Sometimes, PDBFixer will point out some HETATMS that were originally labeled as ATOM.')

    def _defineBoxParams(self, group):
        group.addParam('boxType', params.EnumParam,
                       choices=['Cubic', 'Orthorhombic'],
                       label="Shape of the box: ", default=1,
                       help='Whether to use a Cubic or a Orthorhombic water box')

        group.addParam('sizeType', params.EnumParam,
                       choices=['Absolute', 'Buffer'], display=params.EnumParam.DISPLAY_HLIST,
                       label="System size type: ", default=1,
                       help='Absolute: absolute size of the box (diameter)\n'
                            'Buffer: distance from the solute to the edge of the box')

        line = group.addLine('Box size (nm):',
                             help='Distances of the bounding box (nm).\nIf cubic: a=b=c\n'
                                  'If Orthorhombic: bc=ac=ab=90º')
        line.addParam('distA', params.FloatParam, condition='sizeType == 0',
                      default=5.0, label='a: ')
        line.addParam('distB', params.FloatParam, condition='boxType == 1 and sizeType == 0',
                      default=5.0, label='b: ')
        line.addParam('distC', params.FloatParam, condition='boxType == 1 and sizeType == 0',
                      default=5.0, label='c: ')

        line.addParam('padDist', params.FloatParam, condition='sizeType == 1',
                      default=1.0, label='Buffer distance: ')

    def _defineFFParams(self, group):
        group.addParam('mainForceField', params.EnumParam, choices=GROMACS_LIST,
                       default=GROMACS_AMBER14SB,
                       label='Main Force Field: ',
                       help='Force field applied to the system. Force fields are sets of potential functions and '
                            'parametrized interactions that can be used to study physical systems.')
        group.addParam('waterForceField', params.EnumParam,
                       choices=GROMACS_WATERS_LIST, default=GROMACS_TIP3P,
                       label='Water Force Field: ',
                       help='Force field applied to the waters')

    def _defineIonsParams(self, group):
        group.addParam('placeIons', params.EnumParam, default=1,
                       label='Add ions: ', choices=['None', 'Neutralize', 'Add number'],
                       help='Whether to add ions to the system.'
                            'https://manual.gromacs.org/documentation/2021.5/onlinehelp/gmx-genion.html')

        line = group.addLine('Cation type:', condition='placeIons!=0',
                             help='Type of the cations to add')
        line.addParam('cationType', params.EnumParam, condition='placeIons!=0',
                      label='Cation to add: ', choices=self._cations, default=7,
                      help='Which anion to add in the system')
        line.addParam('cationNum', params.IntParam, condition='placeIons==2',
                      label='Number of cations to add: ',
                      help='Number of cations to add')

        line = group.addLine('Anion type:', condition='placeIons!=0',
                             help='Type of the anions to add')
        line.addParam('anionType', params.EnumParam, condition='placeIons!=0',
                      label='Anions to add: ', choices=self._anions, default=1,
                      help='Which anion to add in the system')
        line.addParam('anionNum', params.IntParam, condition='placeIons==2',
                      label='Number of anions to add: ',
                      help='Number of anions to add')

        group.addParam('addSalt', params.BooleanParam, default=False,
                       condition='placeIons==1',
                       label='Add more salt into the system: ',
                       help='Add more salt into the system')
        group.addParam('saltConc', params.FloatParam, condition='addSalt and placeIons==1',
                       default=0.15,
                       label='Salt concentration (M): ',
                       help='Salt concentration')

    def _defineSSBondsParams(self, group):
        group.addParam('handleSSBonds', params.EnumParam,
                      choices=SSBONDS_OPTIONS,  display=params.EnumParam.DISPLAY_HLIST,
                      default=1, label='Define SS bonds: ',
                      help='How to handle disulfide bonds detected by GROMACS:\n'
                           '*None*: Do not define disulfide bonds\n'
                           '*Automatic*: Define all detected disulfide bonds automatically (SG atoms within 2.0 ± 0.2 Å)\n'
                           '*Manual*: Select which detected bonds to form')

        group.addParam('selectSSBonds', params.StringParam,
                      condition='handleSSBonds==2',
                      label='Selected SS bonds:',
                      help='Comma-separated list of SS bond indices to form (set via wizard)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        if self.usesPreprocessing():
          self._insertFunctionStep(self.preprocessStructureStep)
        if self.inputFrom.get() == LIGAND:
          self._insertFunctionStep(self.parametrizeLigandStep)
        self._insertFunctionStep(self.PDB2GMXStep)
        self._insertFunctionStep(self.editConfStep)
        self._insertFunctionStep(self.solvateStep)
        if self.placeIons.get() != 0:
            self._insertFunctionStep(self.addIonsStep)
        self._insertFunctionStep(self.createOutputStep)

    def preprocessStructureStep(self):
        """PDBFixer repair and/or cleaning of the receptor, written to the file every later step uses"""
        outFile = self.getPreprocessedFile()
        inFile = self.getOriginalReceptorFile()
        if self.usePDBFixer.get():
          self.runPDBFixer(outFile)
          inFile = outFile
        if self.cleansStructure():
          cleanPDB(inFile, outFile, self.waters.get(), self.HETATM.get(), self.getChainsToKeep(),
                   self.getHetResiduesToKeep(), singleModel=0)

    def runPDBFixer(self, outFile):
        addResStr = ' --add-residues' if self.addRes else ''
        repNStdStr = ' --replace-nonstandard' if self.repNonStd else ''
        addAtomsStr = self.getEnumText("addAtoms").lower()

        # PDBFixer reads PDB/mmCIF itself; converting them first would drop the SEQRES records it needs
        inFile = self.getRawInputFile()
        if os.path.splitext(inFile)[1].lower() not in ['.pdb', '.cif', '.mmcif']:
          inFile = self.getOriginalReceptorFile()
        args = f'{inFile} --add-atoms={addAtomsStr}{addResStr}{repNStdStr} --output {outFile}'
        pwchemPlugin.runOPENBABEL(self, 'pdbfixer', args=args, cwd=self._getExtraPath())

        if self.extraClean.get():
          removeStartWithLines(outFile, 'HETATM')

    def parametrizeLigandStep(self):
      mol = self.getSpecifiedMol()
      molFile = os.path.abspath(mol.getPoseFile())
      molFile = self.addHydrogens(molFile)

      kwargs = self.getParameters()
      kwargs['molName'] = mol.getMolName()

      args = f'-i {molFile} -b {kwargs["molName"]} -c {kwargs["chargeMethod"]} ' \
             f'-m {kwargs["multip"]} -a {kwargs["atomType"]} -q {kwargs["qprog"]} -o gmx'
      args += f' -n {kwargs.get("netCharge", getMol2FormalCharge(molFile))}'
      pwchemPlugin.runACPYPE(self, args=args, cwd=self._getExtraPath())


    def addLigandTopo(self, topFile):
      molName = self.getLigandName()

      inStr = '\/forcefield.itp"\n'
      outStr = f'{inStr}\n; Include ligand topology\n#include "{molName}_GMX.itp"\n'
      replaceInFile(topFile, inStr, outStr)

      inStr = '; Compound        #mols\nProtein_chain_A     1'
      outStr = f'{inStr}\n{molName} 1'
      replaceInFile(topFile, inStr, outStr)

    def parseGROFile(self, groFile):
      groDic = {}
      with open(groFile) as f:
        for i, line in enumerate(f):
          if line.strip():
            if i == 0:
              groDic['header'] = line
            elif i == 1:
              groDic['nAtoms'] = line
            else:
              if 'coords' not in groDic:
                groDic['coords'] = ''
              groDic['coords'] += line

      groDic['tail'] = line
      coordsStr = groDic['coords'].replace(groDic['tail'], '')
      groDic['coords'] = coordsStr
      return groDic

    def addLigandCoords(self, recFile, ligFile):
      recDic, ligDic = self.parseGROFile(recFile), self.parseGROFile(ligFile)
      nRec, nLig = int(recDic['nAtoms'].strip()), int(ligDic['nAtoms'].strip())

      with open(recFile, 'w') as f:
        f.write(f"{recDic['header']} {nRec+nLig}\n{recDic['coords']}{ligDic['coords']}{recDic['tail']}")

    def PDB2GMXStep(self):
        inputStructure = self.getInputReceptorFile()
        systemBasename = self.getSystemName()

        addCapsMode = self.getEnumText('addCaps')
        if addCapsMode in GAPS_OPTIONS[1:]:
            mode = 'gaps' if addCapsMode == GAPS_OPTIONS[1] else 'all'
            cappedPdb = os.path.abspath(os.path.join(self._getExtraPath(), f'{systemBasename}_capped.pdb'))
            self.runPymol(self.addCapsPml(inputStructure, cappedPdb, mode), self._getExtraPath())
            self.fixPdbTER(cappedPdb)
            inputStructure = cappedPdb

        # Chains that share residue numbers (e.g. two chains both numbered 1-99) would break the output
        # analysis. Detect that overlap and, only then, renumber the residues continuously.
        if self.chainsHaveOverlappingResNumbers(inputStructure):
            renumberedPdb = os.path.abspath(self._getExtraPath(f'{systemBasename}_renumbered.pdb'))
            inputStructure = self.renumberResiduesPDB(inputStructure, renumberedPdb)

        Waterff = GROMACS_WATERFF_NAME[self.waterForceField.get()]
        Mainff = GROMACS_MAINFF_NAME[self.mainForceField.get()]
        params = (f' pdb2gmx -f {inputStructure} -o {systemBasename}_processed.gro '
                  f'-water {Waterff} -ff {Mainff} -merge all ')

        # Handle disulfide bonds
        ssMode = self.getEnumText('handleSSBonds')
        printfValues = None

        if ssMode != 'None':
            params += ' -ss '
            numBonds = self.countSSBonds(inputStructure, Waterff, Mainff)

            if ssMode == 'Automatic':
                printfValues = ['y'] * int(numBonds)
                self._log.info(f"Automatically accepting {numBonds} disulfide bond(s)")

            elif ssMode == 'Manual':
                selected = {int(x) for x in self.selectSSBonds.get().split(',')} if self.selectSSBonds.get() else set()
                printfValues = ['y' if i in selected else 'n' for i in range(numBonds)]

        try:
            self._runPdb2gmx(params, printfValues=printfValues)

        except:
            self._log.warning('Conversion to gro failed, trying with -ignh flag')
            if os.path.exists(self._getPath(TOPOL_TOP)):
                os.remove(self._getPath(TOPOL_TOP))
            self._runPdb2gmx(params + ' -ignh', printfValues=printfValues)

        if self.inputFrom.get() == LIGAND:
          molName = self.getLigandName()
          groFile, topFile = self._getPath(f'{systemBasename}_processed.gro'), self._getPath(TOPOL_TOP)
          self.addLigandTopo(topFile)
          self.addLigandCoords(groFile, self.getLigandPath(f'{molName}_GMX.gro'))
          shutil.copy(self.getLigandPath(f'{molName}_GMX.itp'), self._getPath(f'{molName}_GMX.itp'))

    def editConfStep(self):
        systemBasename = self.getSystemName()

        boxType = self.getEnumText('boxType').lower() if self.boxType.get() != 1 else 'triclinic'
        params = ' editconf -f %s_processed.gro ' \
                 '-o %s_newbox.gro ' \
                 '-c -bt %s' % (systemBasename, systemBasename, boxType)

        params += self.getDistanceArgs()

        gromacsPlugin.runGromacs(self, 'gmx', params, cwd=self._getPath())

    def solvateStep(self):
        systemBasename = self.getSystemName()

        waterModel = self.getEnumText('waterForceField')
        waterBox = WATER_BOXES.get(waterModel, waterModel)

        params_solvate = ' solvate -cp %s_newbox.gro -cs %s.gro -o %s_solv.gro' \
                         ' -p topol.top' % (systemBasename, waterBox, systemBasename)

        gromacsPlugin.runGromacs(self, 'gmx', params_solvate, cwd=self._getPath())

    def addIonsStep(self):
        systemBasename = self.getSystemName()
        ions_mdp = os.path.abspath(self.buildIonsMDP())

        params_grompp = 'grompp -f %s -c %s_solv.gro -p ' \
                        'topol.top -o ions.tpr' % (ions_mdp, systemBasename)
        if 'gromos' in self.getEnumText('mainForceField'):
            params_grompp += ' -maxwarn 1'
        gromacsPlugin.runGromacsPrintf(self, printfValues=['SOL'],
                                       args=params_grompp, cwd=self._getPath())

        cation, cc = self.parseIon(self.getEnumText('cationType'))
        anion, ac = self.parseIon(self.getEnumText('anionType'))

        genStr = 'genion -s ions.tpr -o %s_solv_ions.gro -p topol.top ' \
                 '-pname %s -nname %s' % (systemBasename, cation, anion)
        if cc == 2:
            genStr += ' -pq {}'.format(cc)
        if ac == 2:
            genStr += ' -nq {}'.format(ac)

        if self.placeIons.get() == 1:
          genStr += ' -neutral '
        elif self.placeIons.get() == 2:
          genStr += ' -np {} -nn {}'.format(self.cationNum.get(), self.anionNum.get())

        if self.addSalt:
          genStr += ' -conc {}'.format(self.saltConc.get())

        gromacsPlugin.runGromacsPrintf(self, printfValues=['SOL'],
                                       args=genStr, cwd=self._getPath())

    def createOutputStep(self):
        systemBasename = self.getSystemName()

        if self.placeIons.get() != 0:
            groBaseName = '%s_solv_ions.gro' % (systemBasename)
        else:
            groBaseName = '%s_solv.gro' % (systemBasename)

        topoPath, groPath, posrePath = self._getPath(TOPOL_TOP), self._getPath(groBaseName), \
                                       self._getPath('posre.itp')

        chainNames = ','.join(self.getModelChains())
        chains, lengthsDic = self.getModelChainsAndLengths()
        lengths = ','.join(str(value) for value in lengthsDic.values())
        
        groSystem = grobj.GromacsSystem(filename=groPath, topoFile=topoPath,
                                        restrFile=posrePath, chainNames=chainNames, chainLengths=lengths,
                                        ff=self.getEnumText('mainForceField'), wff=self.getEnumText('waterForceField'))

        if self.inputFrom.get() == LIGAND:
            molName = self.getLigandName()
            ligTopFile = self._getPath(f'{molName}_GMX.itp')

            # The residue name is read from the ligand topology
            ligName = gromacsPlugin.getLigandResname(ligTopFile)

            groSystem.setLigandID(ligName)
            groSystem.setLigTopologyFile(ligTopFile)
        else:
            ligName = None

        indexFile = gromacsPlugin.firstIndexCreation(self, groSystem, ligandName=ligName, modelChains=chains, chainLengths=lengthsDic)

        groSystem.setIndexFile(indexFile)
        self._defineOutputs(outputSystem=groSystem)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
        mainFF, waterFF = self.getEnumText('mainForceField'), self.getEnumText('waterForceField')

        # Not every force field ships every water model (amber14sb/amber19sb have opc but no tip4p/tip5p,
        # the classic amber ports the other way round), and pdb2gmx -water rejects the ones it lacks.
        supportedWaters = gromacsPlugin.getForceFieldWaterModels(mainFF)
        if supportedWaters and waterFF not in supportedWaters:
            vals.append('Force field {} does not provide the {} water model. It supports: {}.'.format(
              mainFF, waterFF, ', '.join(w for w in supportedWaters if w in GROMACS_WATERS_LIST)))

        try:
            ions = self.getInputIonResidues()
        except Exception:
            ions = set()
        ffResidues = gromacsPlugin.getForceFieldResidues(mainFF)
        missingIons = sorted(ions - ffResidues)
        if missingIons:
            alternatives = [ff for ff in gromacsPlugin.getForceFieldsWithResidues(missingIons)
                            if ff in GROMACS_LIST]
            vals.append('Force field {} has no topology entry for the ion(s) {} present in the input '
                        'structure, so pdb2gmx will not be able to parameterize them.\n{}'.format(
                          mainFF, ', '.join(missingIons),
                          'Installed force fields defining all of them: {}'.format(', '.join(alternatives))
                          if alternatives else 'No installed force field defines all of them: remove them '
                                               'from the input structure or add a custom topology entry.'))

        if not (self.usePDBFixer.get() and self.repNonStd.get()):
            try:
                hetRes = self.getInputHetResidues()
            except Exception:
                hetRes = {}
            missingHet = sorted(resName for resName, nAtoms in hetRes.items()
                                if nAtoms > 1 and resName not in ffResidues)
            if missingHet:
                vals.append('Force field {} has no topology entry for the hetero residue(s) {} kept in the input '
                            'structure, so pdb2gmx will fail on them.\nRemove them with "Remove ligands HETATM" '
                            '(a ligand must be simulated as the docked input ligand) or replace non-standard '
                            'residues with PDBFixer.'.format(mainFF, ', '.join(missingHet)))

        if self.placeIons.get() != 0:
            ionsDic = {'amber': [CA, CL, CS, K, LI, MG, NA, RB, ZN],
                       'gromos': [CA, CL, CU, CU2, MG, NA, ZN],
                       'oplsaa': [BR, CA, CL, CS, F, I, K, LI, NA, RB],
                       'charmm27': [CA, CL, CS, K, MG, NA, ZN],
                       'charmm36': [CA, CL, CS, K, LI, MG, NA, ZN]}
            for key in ionsDic:
                if self.getEnumText('mainForceField').startswith(key):
                    if not self.getEnumText('cationType') in ionsDic[key]:
                        vals.append('{} cation not available for force field {}.\nAvailable ions: {}'.format(
                          self.getEnumText('cationType'), self.getEnumText('mainForceField'), ', '.join(ionsDic[key])
                        ))

                    if not self.getEnumText('anionType') in ionsDic[key]:
                        vals.append('{} anion not available for force field {}.\nAvailable ions: {}'.format(
                          self.getEnumText('anionType'), self.getEnumText('mainForceField'), ', '.join(ionsDic[key])
                        ))
            if self.getEnumText('mainForceField').startswith('gromos') and \
                    self.getEnumText('waterForceField').startswith('tip'):
                vals.append('GROMOS force fields were parametrized for use with SPC water model.'
                            'They will not behave well with TIP models')
        
        return vals

    def _warnings(self):
        warns = []
        mainFF, waterFF = self.getEnumText('mainForceField'), self.getEnumText('waterForceField')
        try:
            ions = self.getInputIonResidues()
        except Exception:
            ions = set()

        mismatched = {}
        for ion in sorted(ions):
            fittedWater = gromacsPlugin.getIonFittedWaterModel(mainFF, ion)
            if fittedWater and fittedWater != waterFF:
                mismatched[ion] = fittedWater

        if mismatched:
            warns.append('{} parameterizes {} for the {} water model, but the system will be solvated with {}, '
                         'so those ions keep hydration parameters fitted for a different water.\n{}'.format(
                           mainFF, ', '.join(f'{ion} ({water})' for ion, water in mismatched.items()),
                           ' / '.join(sorted(set(mismatched.values()))), waterFF,
                           self.getIonWaterRemedy(ions, sorted(set(mismatched.values())), waterFF)))
        return warns

    def getIonWaterRemedy(self, ions, fittedWaters, waterFF):
        """How to make ion parameters and water model agree: switch water when the force field offers the one
        its ions were fitted for, otherwise switch to a force field that fits its ions for this water."""
        mainFF = self.getEnumText('mainForceField')
        ffWaters = gromacsPlugin.getForceFieldWaterModels(mainFF)
        usableWaters = [w for w in fittedWaters if w in GROMACS_WATERS_LIST and w in ffWaters]
        if len(usableWaters) == 1:
            return 'Solvate with {} water instead.'.format(usableWaters[0])

        alternatives = [ff for ff in gromacsPlugin.getForceFieldsForIons(ions, waterFF) if ff in GROMACS_LIST]
        if alternatives:
            return 'Force fields parameterizing these ions for {} water: {}.'.format(waterFF, ', '.join(alternatives))
        return 'No installed force field parameterizes these ions for {} water.'.format(waterFF)

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            if self.usesPreprocessing():
                summary.append(self.getPreprocessSummary())
            summary.append("This protocol has created a processed gro file with Main force field: *%s* and " \
                           "Water Force Field *%s*." % (GROMACS_MAINFF_NAME[self.mainForceField.get()],
                                                      GROMACS_WATERFF_NAME[self.waterForceField.get()]))

        else:
            summary.append("The protocol has not finished.")
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append("This protocol takes a clean pdb file and it uses the "
                           "GROMACS software in order to transform the file into a gromacs format while applying to it "
                           'the force fields for the system and the water molecules. To do so, it calls "gmx pdb2gmx".'
                           'It produces a position restraint file, a topology file and a post-processed structure '
                           'gromacs file.\nThe the protocol runs "gmx editconf" in order to introduce a box from the '
                           'determined size and it runs "gmx solvate" to put water molecules in the box.\nThen the '
                           'command "gmx grompp" is needed to create a tpr file for the command "gmx genion" which will'
                           ' add the selected ion pairs to neutralize the system.\nThen the protocol runs '
                           '"gmx genrestr" in order to create alternative position restriction files which will be used'
                           ' in further steps (position restriction files less restrictive).\nFinally, the output files'
                           ' are created.' )

        return methods

    def _runPdb2gmx(self, params, printfValues=None):
        if printfValues:
            gromacsPlugin.runGromacsPrintf(self, printfValues=printfValues, args=params, cwd=self._getPath())
        else:
            gromacsPlugin.runGromacs(self, 'gmx', params, cwd=self._getPath())

    def getSpecifiedMol(self):
      for mol in self.inputSetOfMols.get():
        if mol.__str__() == self.inputLigand.get():
          return mol.clone()

      raise ValueError(f'Ligand "{self.inputLigand.get()}" is not in the input set of molecules. '
                       f'Available: {self.getInputMolNames()}')

    def getInputMolNames(self):
      mols = self.inputSetOfMols.get()
      return [mol.__str__() for mol in mols] if mols is not None else []

    def getInputReceptorFile(self):
      '''Receptor file the system is built from: the preprocessed one once that step has run'''
      if self.usesPreprocessing() and os.path.exists(self.getPreprocessedFile()):
        return self.getPreprocessedFile()
      return self.getOriginalReceptorFile()

    def usesPreprocessing(self):
      return self.usePDBFixer.get() or self.cleansStructure()

    def cleansStructure(self):
      return self.waters.get() or self.HETATM.get() or self.rchains.get()

    def getChainsToKeep(self):
      '''Model-chain ids selected with the chain wizard (None keeps all), parsed as pwchem's target preparation'''
      if not self.rchains.get():
        return None
      chainJson = json.loads(self.chain_name.get())
      if 'chain' in chainJson:
        return [f'{chainJson["model"]}-{chainJson["chain"].upper().strip()}']
      return [mChain.strip() for mChain in chainJson["model-chain"].upper().split(',')]

    def getHetResiduesToKeep(self):
      return [het.strip() for het in self.het2keep.get().split(',') if het.strip()]

    def getPreprocessSummary(self):
      summary = []
      if self.usePDBFixer.get():
        summary.append('The input structure was repaired with PDBFixer (missing atoms: {}, missing residues: {}, '
                       'non-standard residues replaced: {}).'.format(self.getEnumText('addAtoms'),
                                                                     'yes' if self.addRes else 'no',
                                                                     'yes' if self.repNonStd else 'no'))
      removed = [name for name, flag in [('waters', self.waters.get()), ('HETATM', self.HETATM.get())] if flag]
      if removed:
        kept = self.getHetResiduesToKeep() if self.HETATM.get() else []
        summary.append('Removed from the input structure: {}{}.'.format(
          ', '.join(removed), ' (kept: {})'.format(', '.join(kept)) if kept else ''))
      if self.rchains.get():
        summary.append('Chains kept: {}.'.format(', '.join(self.getChainsToKeep())))
      return ' '.join(summary)

    def getPreprocessedFile(self):
      inName = os.path.splitext(os.path.basename(self.getRawInputFile()))[0]
      return os.path.abspath(self._getExtraPath(inName + '.pdb'))

    def getOriginalReceptorFile(self):
      if self.inputFrom.get() == LIGAND:
        inputStructure = os.path.abspath(self.inputSetOfMols.get().getProteinFile())
      else:
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())

      if not inputStructure.endswith('.pdb'):
        inputPdb = self.getInputPDBFile(inputStructure)
        if not os.path.exists(inputPdb):
          self.convertReceptor2PDB(inputStructure, inputPdb)
        inputStructure = inputPdb
      return inputStructure

    def convertReceptor2PDB(self, proteinFile, oFile):
        _, inExt = os.path.splitext(os.path.basename(proteinFile))
        tmpOFile = '{}.{}.tmp'.format(oFile, uuid.uuid4().hex)
        args = ' -i {} {} -opdb -O {} -d'.format(inExt[1:], os.path.abspath(proteinFile), tmpOFile)
        runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
        os.replace(tmpOFile, oFile)
        return oFile

    def getInputPDBFile(self, proteinFile):
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))
      return os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

    def getSystemName(self):
      return os.path.basename(self.getInputReceptorFile().split(".")[0])

    def addHydrogens(self, inpFile):
      sysbaseName = os.path.basename(inpFile).split('.')[0]
      inpMol2File = os.path.abspath(self._getExtraPath(sysbaseName + '.mol2'))

      if self.isMeekoPdbqt(inpFile):
        # pdbqt has no bond orders: Meeko rebuilds the whole molecule from the SMILES it wrote in the file
        sdfFile = os.path.abspath(self._getTmpPath(sysbaseName + '.sdf'))
        pwchemPlugin.runCondaCommand(self, f'{os.path.abspath(inpFile)} -s {sdfFile}', RDKIT_DIC, 'mk_export.py',
                                     cwd=self._getTmpPath())
        inpFile = sdfFile

      args = f'{os.path.abspath(inpFile)} -h -O {inpMol2File}'
      runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

      replaceInFile(inpMol2File, 'UNL1', 'LIG')
      return inpMol2File

    @staticmethod
    def isMeekoPdbqt(molFile):
      if not molFile.endswith('.pdbqt'):
        return False
      with open(molFile) as f:
        return any(line.startswith('REMARK SMILES') for line in f)

    def getLigandName(self):
      return self.getSpecifiedMol().getMolName()

    def getLigandPath(self, path=''):
      molName = self.getLigandName()
      return self._getExtraPath(f"{molName}.acpype", path)

    def buildIonsMDP(self):
        outFile = self._getPath('ions.mdp')
        ions_mdp_file = "integrator = steep \n" \
                        "emtol = 1000.0 \n" \
                        "emstep = 0.01 \n" \
                        "nsteps = 50000 \n\n" \
                        "nstlist = 10 \n" \
                        "cutoff-scheme = Verlet \n" \
                        "coulombtype = cutoff \n" \
                        "rcoulomb = 1.0 \n" \
                        "rvdw = 1.0 \n" \
                        "pbc = xyz"
        with open(outFile, 'w') as f:
            f.write(ions_mdp_file)
        return outFile

    def parseIon(self, ion):
        if ion[-2].isdigit():
            name, charge = ion[:-2], int(ion[-2])
        else:
            name, charge = ion[:-1], 1
            if name == 'CU':
                name = 'CU1'
        return name, charge

    def getDistanceArgs(self):
        if self.sizeType.get() == 1:
            distArg = ' -d {}'.format(self.padDist.get())
        else:
            if self.boxType.get() == 1:
                distArg = ' -box {} {} {}'.format(self.distA.get(), self.distB.get(), self.distC.get())
            else:
                distArg = ' -box {}'.format(self.distA.get())
        return distArg

    def getModelChains(self):
        inputStructure = self.getInputReceptorFile()
        if not inputStructure.endswith('.pdb'):
            inputPdb = self.getInputPDBFile(inputStructure)
            inputStructure = self.convertReceptor2PDB(inputStructure, inputPdb)

        structureHandler = AtomicStructHandler()
        structureHandler.read(inputStructure)
        structureHandler.getStructure()
        self.dropNonPolymerChains(structureHandler)
        chains, _ = structureHandler.getModelsChains()
        return list(chains[0].keys())

    def getModelChainsAndLengths(self):
        if self.getEnumText('addCaps') in  GAPS_OPTIONS[1:]:
            systemBasename = self.getSystemName()
            inputStructure = os.path.abspath(os.path.join(self._getExtraPath(), f'{systemBasename}_capped.pdb'))
        else:
            inputStructure = self.getInputReceptorFile()

        if not inputStructure.endswith('.pdb'):
            inputPdb = self.getInputPDBFile(inputStructure)
            inputStructure = self.convertReceptor2PDB(inputStructure, inputPdb)

        structureHandler = AtomicStructHandler()
        structureHandler.read(inputStructure)
        structureHandler.getStructure()
        self.dropNonPolymerChains(structureHandler)

        chains, _ = structureHandler.getModelsChains()

        modelChains = list(chains[0].keys())
        # This dictionary will store the count of residues per chain
        lengths = dict(chains[0])

        return modelChains, lengths

    def dropNonPolymerChains(self, structureHandler):
        """Remove the chains holding no protein / nucleic residue, typically a metal ion that the PDB
        converter placed in a chain of its own"""
        for model in structureHandler.structure:
            nonPolymerIds = [chain.id for chain in model if not any(self.isPolymerResidue(res) for res in chain)]
            for chainId in nonPolymerIds:
                model.detach_child(chainId)

    @staticmethod
    def isPolymerResidue(residue):
        resName = residue.get_resname().strip()
        return PDB.is_aa(resName, standard=True) or resName in NUCLEIC_RESNAMES

    def getInputIonResidues(self):
        """Residue names of the free ions and metals left after the preprocessing."""
        return {resName for resName, nAtoms in self.getInputHetResidues().items() if nAtoms == 1}

    def getInputHetResidues(self):
        """Hetero residues (no waters) that reach pdb2gmx, with their number of atoms"""
        if self.usePDBFixer.get() and self.extraClean.get():
            return {}
        inputFile = self.getRawInputFile()
        ext = os.path.splitext(inputFile)[1].lower() if inputFile else ''
        if ext in ['.cif', '.mmcif']:
            parser = PDB.MMCIFParser(QUIET=True)
        elif ext == '.pdb':
            parser = PDB.PDBParser(QUIET=True)
        else:
            return {}

        structure = parser.get_structure('structure', inputFile)
        hetRes = {res.get_resname().strip(): len(list(res.get_atoms()))
                  for res in next(iter(structure)).get_residues() if res.id[0] not in [' ', 'W']}
        if self.HETATM.get():
            hetRes = {resName: nAtoms for resName, nAtoms in hetRes.items()
                      if resName in self.getHetResiduesToKeep()}
        return hetRes

    def getRawInputFile(self):
        """Input structure as given, with no format conversion (usable before the run starts)."""
        inputObj = self.inputSetOfMols if self.inputFrom.get() == LIGAND else self.inputStructure
        if not inputObj.get():
            return None

        inputFile = self.inputSetOfMols.get().getProteinFile() if self.inputFrom.get() == LIGAND \
            else self.inputStructure.get().getFileName()
        return os.path.abspath(inputFile) if inputFile else None

    @staticmethod
    def readGromacsFatalError(logFile):
        """Extract the 'Fatal error' block Gromacs wrote into logFile, so a failing command reports its
        own message instead of just a non-zero exit status."""
        if not os.path.exists(logFile):
            return ''

        with open(logFile) as f:
            lines = f.read().split('\n')

        for i, line in enumerate(lines):
            if line.startswith('Fatal error'):
                return '\n'.join(lines[i:i + 10]).strip()
        return ''

    def runPymol(self, pymolScript, workinDir):
        # run in the background
        self._log.info('Launching: ' + self._getPymol() + pymolScript)
        self.runJob(f'{self._getPymol()} -cq', pymolScript, cwd=workinDir)

    def _getPymol(self):
        return pwchemPlugin.getEnvPath(OPENBABEL_DIC, 'bin/pymol')

    def identifyTermini(self, inputPdb):
        """
        Parses a PDB file to identify chain N/C protein termini and
        internal gaps that require capping.
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", inputPdb)

        result = {'protein_termini': [], 'gaps': []}

        # Collapsing model and chain loops down to one level drops complexity significantly
        for chain in structure.get_chains():
            residues = [r for r in chain if PDB.is_aa(r)]
            if not residues:
                continue

            segments = []
            currentSeqStart = residues[0].id[1]

            for curr, nextRes in zip(residues, residues[1:]):
                resCurr = curr.id[1]
                resNext = nextRes.id[1]

                # Check for a jump in residue numbering (a gap)
                if resNext != resCurr + 1:
                    segments.append({'n': currentSeqStart, 'c': resCurr})
                    currentSeqStart = resNext

            segments.append({'n': currentSeqStart, 'c': residues[-1].id[1]})

            result['protein_termini'].append({
                'chain': chain.id,
                'n_term': segments[0]['n'],
                'c_term': segments[-1]['c']
            })

            for segCurr, segNext in zip(segments, segments[1:]):
                result['gaps'].append({
                    'chain': chain.id,
                    'c_term': segCurr['c'],  # Needs NME
                    'n_term': segNext['n']  # Needs ACE
                })

        return result

    def addCapPmlCommand(self, chain, resi, atom, capType):
        return [
            f"select tmp_target, /protein//{chain}/{resi}/{atom}",
            "edit tmp_target",
            f"/editor.attach_amino_acid('pk1', '{capType}')"
        ]

    def addCapsPml(self, inputPdb, outputPdb, mode='gaps'):
        data = self.identifyTermini(inputPdb)

        pmlLines = [
            "reinitialize",
            f"load {inputPdb}, protein",
            "remove hydro",
            "hide all",
            "show sticks, protein"
        ]
        if 'gaps' in data and data['gaps']:
            for gap in data['gaps']:
                # Adds NME on the C-term and ACE on the N-term of gaps
                pmlLines.append(self.removeOXTCommand(gap['chain'], gap['c_term']))
                pmlLines.extend(self.addCapPmlCommand(gap['chain'], gap['c_term'], 'C', 'nme'))
                pmlLines.extend(self.addCapPmlCommand(gap['chain'], gap['n_term'], 'N', 'ace'))

        if mode == 'all':
            for term in data['protein_termini']:
                # Adds NME on the C-term and ACE on the N-term of chain termini
                pmlLines.append(self.removeOXTCommand(term['chain'], term['c_term']))
                pmlLines.extend(self.addCapPmlCommand(term['chain'], term['n_term'], 'N', 'ace'))
                pmlLines.extend(self.addCapPmlCommand(term['chain'], term['c_term'], 'C', 'nme'))

        pmlLines.extend([
            "remove hydro",
            "sort protein",
            f"save {outputPdb}, protein",
            "quit"
        ])

        capScript = os.path.abspath(os.path.join(self._getExtraPath(),"capping_script.pml"))
        with open(capScript, "w") as f:
            f.write("\n".join(pmlLines))

        print(f"PML script for adding caps: {capScript} for mode: {mode}")
        return capScript

    def removeOXTCommand(self, chain, resi):
        return f"remove /protein//{chain}/{resi}/OXT"

    def fixPdbTER(self, pdbPath):
        with open(pdbPath, 'r') as f:
            lines = f.readlines()

        cleanLines = [line for line in lines if not line.startswith("TER")]
        fixedLines = []
        numLines = len(cleanLines)

        for i, line in enumerate(cleanLines):
            fixedLines.append(line)

            # Guard Clause: Skip any line that isn't an NME ATOM record
            if not line.startswith("ATOM") or line[17:20].strip() != "NME":
                continue

            # If it's the very last line, it needs a TER
            if i + 1 >= numLines:
                fixedLines.append("TER\n")
                continue

            # Otherwise, check if the next atom belongs to a different residue
            nextLine = cleanLines[i + 1]
            if nextLine.startswith("ATOM") and nextLine[22:26].strip() != line[22:26].strip():
                fixedLines.append("TER\n")

        with open(pdbPath, 'w') as f:
            f.writelines(fixedLines)

    def chainsHaveOverlappingResNumbers(self, inPdb):
        """Return True if two or more chains share at least one residue number (e.g. two chains both
        numbered 1-99). """
        parser = PDB.PDBParser(QUIET=True)
        model = next(iter(parser.get_structure('protein', inPdb)))
        seen = set()
        for chain in model:
            # res.id is (hetero-flag, residue_number, insertion_code); res.id[1] is the residue number
            chainNums = {res.id[1] for res in chain}
            if seen & chainNums:
                return True
            seen |= chainNums
        return False

    def renumberResiduesPDB(self, inPdb, outPdb):
        """Renumber the residues of a PDB continuously across all chains so that, after merging the chains
        with "gmx pdb2gmx -merge all", the resulting .gro does not contain repeated residue numbers.
        """
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', inPdb)
        # Renumber only the first model (the one pdb2gmx will read)
        model = next(iter(structure))
        residues = list(model.get_residues())
        startNum = residues[0].id[1]
        # Two loops are required so we first park every residue at a temporary number
        # (100000 + i), guaranteed free and far above any real residue number, and then assign the final.
        for i, res in enumerate(residues):
            res.id = (res.id[0], 100000 + i, ' ')
        for offset, res in enumerate(residues):
            res.id = (res.id[0], startNum + offset, ' ')

        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(outPdb)
        self._log.info(f'Residues renumbered continuously into {outPdb}')
        return outPdb

    def countSSBonds(self, inputStructure, waterff='spc', mainff='amber03'):
        """
        Run pdb2gmx with -ss to count how many SS bonds GROMACS detects.
        Returns the count without forming any bonds.
        """
        maxBonds = self.maximumSSbonds(inputStructure)
        outputFile = os.path.abspath(self._getTmpPath('pdb2gmx_output.log'))
        params = f'pdb2gmx -f {inputStructure} -o test.gro -water {waterff} -ff {mainff} -merge all -ss -ignh' \
               f' > {outputFile} 2>&1'

        printfValues =['n'] * maxBonds
        try:
            gromacsPlugin.runGromacsPrintf(self, printfValues, params, cwd=self._getTmpPath())
        except Exception:
            raise RuntimeError('gmx pdb2gmx failed while detecting disulfide bonds.\n{}'.format(
              self.readGromacsFatalError(outputFile)))

        # Count "Link CYS..." lines
        with open(outputFile, 'r') as f:
            output = f.read()
        count = output.count('(y/n) ?')
        return count

    def maximumSSbonds(self, inputStructure):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", inputStructure)
        cysCount = sum(1 for residue in structure.get_residues()
                        if residue.get_resname() in ['CYS', 'CYX'])
        # Maximum possible bonds is = n*(n-1)/2
        maxBonds = (cysCount * (cysCount - 1)) // 2
        return maxBonds

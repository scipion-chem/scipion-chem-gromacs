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
import os, subprocess, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.convert import AtomicStructHandler

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import runOpenBabel
from pwchem.protocols import ProtocolLiganParametrization

from gromacs import Plugin as gromacsPlugin
import gromacs.objects as grobj
from gromacs.constants import *

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

GROMACS_LIST = [GROMACS_MAINFF_NAME[GROMACS_AMBER03], GROMACS_MAINFF_NAME[GROMACS_AMBER94],
                GROMACS_MAINFF_NAME[GROMACS_AMBER96], GROMACS_MAINFF_NAME[GROMACS_AMBER99],
                GROMACS_MAINFF_NAME[GROMACS_AMBER99SB], GROMACS_MAINFF_NAME[GROMACS_AMBERSB_ILDN],
                GROMACS_MAINFF_NAME[GROMACS_AMBERGS],
                GROMACS_MAINFF_NAME[GROMACS_CHARMM27], GROMACS_MAINFF_NAME[GROMACS_CHARMM36],
                GROMACS_MAINFF_NAME[GROMACS_GROMOS43A1], GROMACS_MAINFF_NAME[GROMACS_GROMOS43A2],
                GROMACS_MAINFF_NAME[GROMACS_GROMOS45A3], GROMACS_MAINFF_NAME[GROMACS_GROMOS53A5],
                GROMACS_MAINFF_NAME[GROMACS_GROMOS53A6], GROMACS_MAINFF_NAME[GROMACS_GROMOS54A7],
                GROMACS_MAINFF_NAME[GROMACS_OPLSAA]]

GROMACS_SPC = 0
GROMACS_SPCE = 1
GROMACS_TIP3P = 2
GROMACS_TIP4P = 3
GROMACS_TIP5P = 4

GROMACS_WATERFF_NAME = dict()
GROMACS_WATERFF_NAME[GROMACS_SPC] = 'spc'
GROMACS_WATERFF_NAME[GROMACS_SPCE] = 'spce'
GROMACS_WATERFF_NAME[GROMACS_TIP3P] = 'tip3p'
GROMACS_WATERFF_NAME[GROMACS_TIP4P] = 'tip4p'
GROMACS_WATERFF_NAME[GROMACS_TIP5P] = 'tip5p'

GROMACS_WATERS_LIST = [GROMACS_WATERFF_NAME[GROMACS_SPC], GROMACS_WATERFF_NAME[GROMACS_SPCE],
GROMACS_WATERFF_NAME[GROMACS_TIP3P], GROMACS_WATERFF_NAME[GROMACS_TIP4P],
GROMACS_WATERFF_NAME[GROMACS_TIP5P]]

STRUCTURE, LIGAND = 0, 1

def replaceInFile(file, inStr, repStr):
  inStr, repStr = inStr.replace('\n', '\\n'), repStr.replace('\n', '\\n')
  subprocess.check_call(f"sed -i -z 's/{inStr}/{repStr}/g' {file}", shell=True)
  return file


class GromacsSystemPrep(ProtocolLiganParametrization):
    """
    This protocol will start a Molecular Dynamics preparation. It will create the system
    and the topology, structure, and position restriction files

    It is necessary to insert a cleaned PDB strucrture from Protocol Import Atomic Structure
    or other similar protocols.
    """
    _label = 'system preparation'
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

        form.addParam('inputFrom', params.EnumParam, default=STRUCTURE,
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
        self._defineACPYPEparams(form, condition=f'inputFrom=={LIGAND}')

        group = form.addGroup('Boundary box')
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

        form.addSection('Force Field')
        group = form.addGroup('Force field')
        group.addParam('mainForceField', params.EnumParam, choices=GROMACS_LIST,
                       default=GROMACS_AMBER03,
                       label='Main Force Field: ',
                       help='Force field applied to the system. Force fields are sets of potential functions and '
                            'parametrized interactions that can be used to study physical systems.')
        group.addParam('waterForceField', params.EnumParam,
                       choices=GROMACS_WATERS_LIST, default=GROMACS_TIP3P,
                       label='Water Force Field: ',
                       help='Force field applied to the waters')

        group = form.addGroup('Ions')
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

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        if self.inputFrom.get() == LIGAND:
          self._insertFunctionStep(self.parametrizeLigandStep)
        self._insertFunctionStep(self.PDB2GMXStep)
        self._insertFunctionStep(self.editConfStep)
        self._insertFunctionStep(self.solvateStep)
        if self.placeIons.get() != 0:
            self._insertFunctionStep(self.addIonsStep)
        self._insertFunctionStep(self.createOutputStep)

    def parametrizeLigandStep(self):
      mol = self.getSpecifiedMol()
      molFile = os.path.abspath(mol.getPoseFile())
      molFile = self.addHydrogens(molFile)

      kwargs = self.getParameters()
      kwargs['molName'] = mol.getMolName()

      args = f'-i {molFile} -b {kwargs["molName"]} -c {kwargs["chargeMethod"]} ' \
             f'-m {kwargs["multip"]} -a {kwargs["atomType"]} -q {kwargs["qprog"]} -o gmx'
      if 'netCharge' in kwargs:
        args += f' -n {kwargs["netCharge"]}'
      pwchemPlugin.runACPYPE(self, args=args, cwd=self._getExtraPath())


    def addLigandTopo(self, topFile):
      molName = self.getLigandName()

      inStr = '\/forcefield.itp"\n'
      outStr = f'{inStr}\n; Include ligand topology\n#include "{molName}_GMX.itp"\n'
      replaceInFile(topFile, inStr, outStr)

      emptyStr = ' ' * (20-len(molName))
      inStr = '; Compound        #mols\nProtein_chain_A     1'
      outStr = f'{inStr}\n{molName}{emptyStr}1'
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

      Waterff = GROMACS_WATERFF_NAME[self.waterForceField.get()]
      Mainff = GROMACS_MAINFF_NAME[self.mainForceField.get()]
      params = f' pdb2gmx -f {inputStructure} -o {systemBasename}_processed.gro ' \
               f'-water {Waterff} -ff {Mainff} -merge all'
      # todo: managing several chains (restrictions, topologies...) instead of merging them
      try:
          gromacsPlugin.runGromacs(self, 'gmx', params, cwd=self._getPath())
      except:
          print('Conversion to gro failed, trying to convert it ignoring the current hydrogens')
          os.remove(self._getPath('topol.top'))
          params += ' -ignh'
          gromacsPlugin.runGromacs(self, 'gmx', params, cwd=self._getPath())

      if self.inputFrom.get() == LIGAND:
        molName = self.getLigandName()
        groFile, topFile = self._getPath(f'{systemBasename}_processed.gro'), self._getPath('topol.top')
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
        if waterModel in ['spc', 'spce', 'tip3p']:
            waterModel = 'spc216'

        params_solvate = ' solvate -cp %s_newbox.gro -cs %s.gro -o %s_solv.gro' \
                         ' -p topol.top' % (systemBasename, waterModel, systemBasename)

        gromacsPlugin.runGromacs(self, 'gmx', params_solvate, cwd=self._getPath())

    def addIonsStep(self):
        systemBasename = self.getSystemName()
        ions_mdp = os.path.abspath(self.buildIonsMDP())

        params_grompp = 'grompp -f %s -c %s_solv.gro -p ' \
                        'topol.top -o ions.tpr' % (ions_mdp, systemBasename)
        if 'gromos' in self.getEnumText('mainForceField'):
            params_grompp += ' -maxwarn 1'
        gromacsPlugin.runGromacsPrintf(printfValues=['SOL'],
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

        gromacsPlugin.runGromacsPrintf(printfValues=['SOL'],
                                       args=genStr, cwd=self._getPath())

    def createOutputStep(self):
        systemBasename = self.getSystemName()

        if self.placeIons.get() != 0:
            groBaseName = '%s_solv_ions.gro' % (systemBasename)
        else:
            groBaseName = '%s_solv.gro' % (systemBasename)

        topoPath, groPath, posrePath = self._getPath('topol.top'), self._getPath(groBaseName), \
                                       self._getPath('posre.itp')

        chainNames = ','.join(self.getModelChains())

        groSystem = grobj.GromacsSystem(filename=groPath, topoFile=topoPath,
                                        restrFile=posrePath, chainNames=chainNames,
                                        ff=self.getEnumText('mainForceField'), wff=self.getEnumText('waterForceField'))
        if self.inputFrom.get() == LIGAND:
          molName = self.getLigandName()
          groSystem.setLigandName(molName)
          groSystem.setLigandTopologyFile(self._getPath(f'{molName}_GMX.itp'))

        self._defineOutputs(outputSystem=groSystem)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        vals = []
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

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
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

    def getSpecifiedMol(self):
      myMol = None
      for mol in self.inputSetOfMols.get():
        if mol.__str__() == self.inputLigand.get():
          myMol = mol.clone()
          break
      if myMol == None:
        print('The input ligand is not found')
        return None
      else:
        return myMol

    def getInputReceptorFile(self):
      if self.inputFrom.get() == LIGAND:
        inputStructure = os.path.abspath(self.inputSetOfMols.get().getProteinFile())
      else:
        inputStructure = os.path.abspath(self.inputStructure.get().getFileName())

      if not inputStructure.endswith('.pdb'):
        inputStructure = self.getInputPDBFile(inputStructure)
        if not os.path.exists(inputStructure):
          inputStructure = self.convertReceptor2PDB(inputStructure)
      return inputStructure

    def getInputPDBFile(self, proteinFile):
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))
      return os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

    def getSystemName(self):
      return os.path.basename(self.getInputReceptorFile().split(".")[0])

    def addHydrogens(self, inpFile):
      sysbaseName = os.path.basename(inpFile).split('.')[0]
      tmpFile = os.path.abspath(self._getTmpPath(sysbaseName + '.pdb'))
      inpMol2File = os.path.abspath(self._getExtraPath(sysbaseName + '.mol2'))

      args = f'{os.path.abspath(inpFile)} -O {tmpFile}'
      runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

      args = f'{os.path.abspath(tmpFile)} -h -O {inpMol2File}'
      runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

      replaceInFile(inpMol2File, 'UNL1', 'LIG')
      return inpMol2File

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

    def convertReceptor2PDB(self, proteinFile):
        _, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = self.getInputPDBFile(proteinFile)
        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
        runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        return oFile

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
          inputStructure = self.convertReceptor2PDB(inputStructure)

        structureHandler = AtomicStructHandler()
        structureHandler.read(inputStructure)
        structureHandler.getStructure()
        chains, _ = structureHandler.getModelsChains()
        return list(chains[0].keys())

# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Pedro Febrer Martinez (pedrofebrer98@gmail.com)
# *
# * your institution
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
This module will prepare the system for the simumlation
"""
import os
from os.path import abspath, relpath
from pyworkflow.protocol import Protocol, params, Integer
from pyworkflow.utils import Message, runJob, createLink
import pwem.objects as emobj
import gromacs.objects as grobj
from pwem.protocols import EMProtocol
import shutil
import json

ION_NA = 0
ION_K = 1

ION_NAME = dict()
ION_NAME[ION_NA] = 'NA'
ION_NAME[ION_K] = 'K'

ION_LIST = [ION_NAME[ION_NA], ION_NAME[ION_K]]

GROMACS_AMBER03 = 0
GROMACS_AMBER94 = 1
GROMACS_AMBER96 = 2
GROMACS_AMBER99 = 3
GROMACS_AMBER99SB = 4
GROMACS_AMBERSB_ILDN = 5
GROMACS_AMBERGS = 6
GROMACS_CHARMM27 = 7
GROMACS_GROMOS43A1 = 8
GROMACS_GROMOS43A2 = 9
GROMACS_GROMOS45A3 = 10
GROMACS_GROMOS53A5 = 11
GROMACS_GROMOS53A6 = 12
GROMACS_GROMOS54A7 = 13
GROMACS_OPLSAA = 14

GROMACS_MAINFF_NAME = dict()
GROMACS_MAINFF_NAME[GROMACS_AMBER03] = 'amber03'
GROMACS_MAINFF_NAME[GROMACS_AMBER94] = 'amber94'
GROMACS_MAINFF_NAME[GROMACS_AMBER96] = 'amber96'
GROMACS_MAINFF_NAME[GROMACS_AMBER99] = 'amber99'
GROMACS_MAINFF_NAME[GROMACS_AMBER99SB] = 'amber99sb'
GROMACS_MAINFF_NAME[GROMACS_AMBERSB_ILDN] = 'amber99sb-ildn'
GROMACS_MAINFF_NAME[GROMACS_AMBERGS] = 'amberGS'
GROMACS_MAINFF_NAME[GROMACS_CHARMM27] = 'charmm27'
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
GROMACS_MAINFF_NAME[GROMACS_AMBERGS], GROMACS_MAINFF_NAME[GROMACS_CHARMM27],
GROMACS_MAINFF_NAME[GROMACS_GROMOS43A1], GROMACS_MAINFF_NAME[GROMACS_GROMOS43A2],
GROMACS_MAINFF_NAME[GROMACS_GROMOS45A3], GROMACS_MAINFF_NAME[GROMACS_GROMOS53A5],
GROMACS_MAINFF_NAME[GROMACS_GROMOS53A6], GROMACS_MAINFF_NAME[GROMACS_GROMOS54A7],
GROMACS_MAINFF_NAME[GROMACS_OPLSAA]]

GROMACS_SPC = 0
GROMACS_SPCE = 1
GROMACS_TIP3P = 2
GROMACS_TIP4P = 3
GROMACS_TIP5P = 4
GROMACS_NONE = 5

GROMACS_WATERFF_NAME = dict()
GROMACS_WATERFF_NAME[GROMACS_SPC] = 'spc'
GROMACS_WATERFF_NAME[GROMACS_SPCE] = 'spce'
GROMACS_WATERFF_NAME[GROMACS_TIP3P] = 'tip3p'
GROMACS_WATERFF_NAME[GROMACS_TIP4P] = 'tip4p'
GROMACS_WATERFF_NAME[GROMACS_TIP5P] = 'tip5p'
GROMACS_WATERFF_NAME[GROMACS_NONE] = 'none'

GROMACS_WATERS_LIST = [GROMACS_WATERFF_NAME[GROMACS_SPC], GROMACS_WATERFF_NAME[GROMACS_SPCE],
GROMACS_WATERFF_NAME[GROMACS_TIP3P], GROMACS_WATERFF_NAME[GROMACS_TIP4P],
GROMACS_WATERFF_NAME[GROMACS_TIP5P],
GROMACS_WATERFF_NAME[GROMACS_NONE]]

class GromacsSystemPrep(EMProtocol):
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

    # -------------------------- DEFINE constants ----------------------------


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('UsePDBFile', params.PointerParam,
                      label=" PDB to use", allowsNull=False,
                      important=True,
                      pointerClass='AtomStruct',
                      help='This PDB file will be used to do pdb2gmx')

        form.addParam('mainForceField', params.EnumParam,
                      choices=GROMACS_LIST,
                      default=GROMACS_AMBER03,
                      label='Main Force Field', important=True,
                      help='Force field applied to the system. Force fields are sets of potential functions and '
                           'parametrized interactions that can be used to study physical systems.')

        form.addParam('waterForceFieldList', params.EnumParam,
                      choices=GROMACS_WATERS_LIST,
                      default=GROMACS_TIP3P,
                      label='Water Force Field', important=True,
                      help='Force field applied to the waters')

        form.addParam('boxsize', params.FloatParam,
                      label="Insert distance from protein to box (nm)",
                      default=1.0, important=True,
                      help='This value should be between 1.0 and 1.5 nm. The higher the value, the higher the '
                           'computational cost. It depends on the force field chosen.'\
                           'The box is a cube.')

        form.addParam('energy', params.IntParam,
                      label="Insert energy restriction",
                      default=1000, important=True,
                      allowsnull=False,
                      help='Force constant for position restraints applied to heavy atoms in the system.' 
                           '\nThis value should be between 1000 and 50000.')

        form.addParam('addIonCharges', params.EnumParam,
                      label="Add Ion types",
                      choices=('NA+/Cl-','K+/Cl-'),
                      allowsNull=False, important=True,
                      default=1,
                      help='Ions which will be added to the system in order to make it neutral. Some random water'
                           'molecules will be replaced by some ions. ')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('getPDB2GMXParams')
        self._insertFunctionStep('getEditconfParams')
        self._insertFunctionStep('getEditconfParams')
        self._insertFunctionStep('getSolvateParams')
        self._insertFunctionStep('getGromppParams')
        self._insertFunctionStep('getGenionParams')
        self._insertFunctionStep('getGenrestrParams')
        self._insertFunctionStep('createOutputStep')

    def getPDB2GMXParams(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])
        print(OutputPDBFile)
        Waterff = GROMACS_WATERFF_NAME[self.waterForceFieldList.get()]
        Mainff = GROMACS_MAINFF_NAME[self.mainForceField.get()]
        energy = self.energy.get()
        program = os.path.join("",'/usr/local/gromacs/bin/gmx')
        print(program)
        params = ' pdb2gmx -f %s ' \
                 '-o %s_processed.gro ' \
                 '-water %s ' \
                 '-ff %s ' \
                 '-posrefc %d' % (UsePDBFile, OutputPDBFile, Waterff, Mainff, energy)
        self.runJob(program, params, cwd=self._getPath())

    def getEditconfParams(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        params = ' editconf -f %s_processed.gro ' \
                 '-o %s_newbox.gro ' \
                 '-c -d %s -bt cubic' % (OutputPDBFile, OutputPDBFile, self.boxsize.get())

        self.runJob(program, params, cwd=self._getPath())

    def getSolvateParams(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')

        params_solvate = ' solvate -cp %s_newbox.gro -cs spc216.gro -o %s_solv.gro' \
                         ' -p topol.top' % (OutputPDBFile, OutputPDBFile)

        self.runJob(program, params_solvate, cwd=self._getPath())

    def getGromppParams(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        working_dir = self.getWorkingDir()
        ions_mdp = '%s/ions.mdp' % (working_dir)
        ions_mdp_file = "integrator = steep \n" \
                        "emtol = 1000.0 \n" \
                        "emstep = 0.01 \n" \
                        "nsteps = 50000 \n" \
                        "\n" \
                        "nstlist = 1 \n" \
                        "cutoff-scheme = Verlet \n" \
                        "ns_type = grid \n" \
                        "coulombtype = cutoff \n" \
                        "rcoulomb = 1.0 \n" \
                        "rvdw = 1.0 \n" \
                        "pbc = xyz"
        f = open(ions_mdp, "w")
        f.write(ions_mdp_file)
        f.close()

        params_grompp = 'grompp -f ions.mdp -c %s_solv.gro -p ' \
                        'topol.top -o ions.tpr' % (OutputPDBFile)
        self.runJob(program, params_grompp, cwd=self._getPath())

    def getGenionParams(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])
        program = os.path.join("", 'printf "13" | /usr/local/gromacs/bin/gmx')
        ION = ION_NAME[self.addIonCharges.get()]

        params_genion = 'genion -s ions.tpr ' \
                 '-o %s_solv_ions.gro ' \
                 '-p topol.top -pname %s -nname CL -neutral' % (OutputPDBFile, ION)
        self.runJob(program, params_genion, cwd=self._getPath())

    def getGenrestrParams(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])
        energy = (self.energy.get())/3.2
        program = os.path.join("", 'printf "2" | /usr/local/gromacs/bin/gmx')
        params_genrestr = 'genrestr -f %s_solv_ions.gro -o posre_low.itp -fc %d %d %d' % (OutputPDBFile,
                                                                                                energy, energy, energy)
        self.runJob(program, params_genrestr, cwd=self._getPath())
        program = "sed "
        sed_params = """-i '/; Include Position restraint file/a #ifdef POSRES_LOW' topol.top"""
        self.runJob(program, sed_params, cwd=self._getPath())
        sed_params = """-i '/#ifdef POSRES_LOW/a #include "posre_low.itp"' topol.top"""
        self.runJob(program, sed_params, cwd=self._getPath())
        sed_params = """-i '/#include "posre_low.itp"/a #endif' topol.top"""
        self.runJob(program, sed_params, cwd=self._getPath())


    def createOutputStep(self):
        UsePDBFile = os.path.abspath(self.UsePDBFile.get().getFileName())
        OutputPDBFile = os.path.basename(UsePDBFile.split(".")[0])

        gro_baseName = '%s_solv_ions.gro' % (OutputPDBFile)
        topol_baseName = 'topol.top'
        posre_baseName = 'posre.itp'
        posre_h_baseName = 'posre_low.itp'

        topol_localPath = relpath(abspath(self._getPath(topol_baseName)))
        gro_localPath = relpath(abspath(self._getPath(gro_baseName)))
        posre_localPath = relpath(abspath(self._getPath(posre_baseName)))
        posre_h_localPath = relpath(abspath(self._getPath(posre_h_baseName)))

        gro_files = grobj.GroFiles(gro=gro_localPath,
                                    topol=topol_localPath,
                                    posre_h=posre_h_localPath,
                                    posre=posre_localPath)

        self._defineOutputs(outputGroFiles=gro_files)
        self._defineSourceRelation(self.UsePDBFile, gro_files)

    # --------------------------- INFO functions -----------------------------------


    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            summary.append("This protocol has created a processed gro file with Main force field: *%s* and " \
                           "Water Force Field *%s*." % (GROMACS_MAINFF_NAME[self.mainForceField.get()],
                                                      GROMACS_WATERFF_NAME[self.waterForceFieldList.get()]))

        else:
            summary.append("The protocol has not finished.")
        return summary
    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append("This protocol takes a clean pdb file and it uses the "
                           "GROMACS software in order to transform the file into a gromacs format while applying to it "
                           'the force fields for the system and the water molecules. To do so, it calls "gmx pdb2gmx".'
                           'It produces a position restrain file, a topology file and a post-processed structure '
                           'gromacs file.\nThe the protocol runs "gmx editconf" in order to introduce a box from the '
                           'determined size and it runs "gmx solvate" to put water molecules in the box.\nThen the '
                           'command "gmx grompp" is needed to create a tpr file for the command "gmx genion" which will'
                           ' add the selected ion pairs to neutralize the system.\nThen the protocol runs '
                           '"gmx genrestr" in order to create alternative position restriction files which will be used'
                           ' in further steps (position restriction files less restrictive).\nFinally, the output files'
                           ' are created.' )

        return methods

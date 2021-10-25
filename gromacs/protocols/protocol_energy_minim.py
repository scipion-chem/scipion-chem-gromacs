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
This module will perform energy minimizations for the system
"""
import os
from os.path import exists, basename, abspath, relpath, join
import pyworkflow.utils as pwutils
from pyworkflow.protocol import Protocol, params, Integer, EnumParam, StringParam, FloatParam, BooleanParam
from pyworkflow.utils import Message, runJob, createLink
import pwem.objects as emobj
import gromacs.objects as grobj
from gromacs.objects import *
from pyworkflow.protocol.params import (LEVEL_ADVANCED, LEVEL_NORMAL, USE_GPU)
from pwem.protocols import EMProtocol

GROMACS_SPC = 0
GROMACS_SPCE = 1
GROMACS_TIP3P = 2


GROMACS_WATERFF_NAME = dict()
GROMACS_WATERFF_NAME[GROMACS_SPC] = 'steep'
GROMACS_WATERFF_NAME[GROMACS_SPCE] = 'mp'
GROMACS_WATERFF_NAME[GROMACS_TIP3P] = 'md-vv'

INTEGRATOR_LIST = [[GROMACS_WATERFF_NAME[GROMACS_SPC], GROMACS_WATERFF_NAME[GROMACS_SPCE],
GROMACS_WATERFF_NAME[GROMACS_TIP3P]]]
class GromacsEnergyMinimization(EMProtocol):
    """
    This protocol will perform energy minimization on the system previosly prepared by the protocol "system prepartion".
    This step is necessary to energy minize the system in order to avoid unwanted conformations.
    """
    _label = 'energy minimization'
    EM_OPTIONS = ['Hydrogen+Water+System', 'Water oxygen+System', 'System']
    EM_PARAMS_OPTIONS = ['steep', 'cg']
    # -------------------------- DEFINE constants ----------------------------


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)

        form.addHidden(params.USE_GPU, params.BooleanParam,
                      label='Use GPU for execution',
                      default=False,
                      help="GPU may have several cores. Set it one if "
                           "you don't know what we are talking about but you have a GPU."
                           "For DARC, first core index is 1, second 2, and so on. Write 0 if you do not want"
                           "to use GPU")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addParam('ScipionOrLocalParams', params.BooleanParam,
                      label="Use default parameters",
                      default=True,
                      expertLevel=LEVEL_NORMAL,
                      condition="expertLevel == LEVEL_NORMAL",
                      help="'Yes' = It uses the default params for energy minimization: \n"
                        "integrator = steep ; Algorithm (steep = steepest descent minimization)\n"
                        "emtol = 1000.0 ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm\n"
                        "emstep = 0.01 ; Minimization step size\n"
                        "nsteps = 50000 ; Maximum number if (minimization steps to perform)\n"
                        "nstlist = 1 ; Frequency to update the neighbor list and long range forces\n"
                        "cutoff-scheme = Verlet ; Buffered neighbor searching\n"
                        "ns_type = grid ; Method to determine neighbor list (simple, grid)\n"
                        "coulombtype = PME ; Treatment of long range electrostatic interactions\n"
                        "rcoulomb = 1.0 ; Shor-range electrostatic cut-off\n"
                        "rvdw = 1.0 ; Shord-range Van der Waals cut-off\n"
                        "pbc = xyz ; Periodic Boundary Conditions in all 3 dimensions\n"
                        "\n"
                      "'No' == It uses the params from an input parameters file.")

        form.addParam('EMSteeps', params.EnumParam,
                      label='Energy minimization steps',
                      expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=self.EM_OPTIONS,
                      default=0,
                      allowsnull=False,
                      help='Energy minimization steps which will take place in the protocol.\n*"Hydrogen+Water+System"*'
                           ' will perform 3 minimization steps: First, it will minimize all hydrogens in the system, '
                           'applying energy restrictions to the rest of the atoms; then it will minimize all oxygens '
                           'from water molecules applying energy restrictions only to non-hydrogen atoms from the '
                           'protein and it finally will minimize the system without any restrictions.\n'
                           '*"Water oxygen+System"* will perform 2 minimization steps: First, it will minimize all water '
                           'oxygens and hydrogens, with restrictions to the heavy atoms of the proteins and finally it '
                           'will minimize the whole system without any restrictions.\n'
                           '*"System"* will mimimize only the whole system without restrictions.\n')

        form.addParam('EMAlgorithmHydrogen', params.EnumParam,
                      label='Algorithm for hydrogen energy minimization',
                      expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_HLIST,
                      condition='EMSteeps == 0',
                      choices=self.EM_PARAMS_OPTIONS,
                      allowsnull=False,
                      default=0,
                      help='Algorithm used to run the energy minimization for the hydrogens.\n*"steep"* will use the '
                           'steepest descent algorithm.\n*"cg"* will use a conjugate gradient algorithm which will also'
                           ' do steepest descent every X steeps.\nSteepest descent is more accurate and time consuming'
                           ' than conjugate gradient. It is recommended to perform steepest descent for the first 2 '
                           'minimizations and conjugate gradient for the last one. In case you only perform 1 '
                           'minimzation it is recommended to perform steepest descent which is more accurate.')

        form.addParam('EMNstcgsteepHydrogen', params.IntParam,
                      label='Nstcgsteep for CG algorithm in Hydrogen minimization:',
                      expertLevel=LEVEL_ADVANCED,
                      condition='EMAlgorithmHydrogen == 1',
                      default=1000,
                      allowsnull=False,
                      help="Every X amount of steps of conjugate gradient algorithm it will be performed and steepest "
                           "descent algorithm.")

        form.addParam('EMAlgorithmWater', params.EnumParam,
                      label='Algorithm for waters energy minimization',
                      expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_HLIST,
                      condition='EMSteeps == 0 or EMSteeps == 1',
                      choices=self.EM_PARAMS_OPTIONS,
                      allowsnull=False,
                      default=0,
                      help='Algorithm used to run the energy minimization for the hydrogens.\n*"steep"* will use the '
                           'steepest descent algorithm.\n*"cg"* will use a conjugate gradient algorithm which will also'
                           ' do steepest descent every X steeps.\nSteepest descent is more accurate and time consuming'
                           ' than conjugate gradient. It is recommended to perform steepest descent for the first 2 '
                           'minimizations and conjugate gradient for the last one. In case you only perform 1 '
                           'minimzation it is recommended to perform steepest descent which is more accurate.')

        form.addParam('EMNstcgsteepWater', params.IntParam,
                      label='Nstcgsteep for CG algorithm in Water minimization:',
                      expertLevel=LEVEL_ADVANCED,
                      condition='EMAlgorithmWater == 1',
                      default=1000,
                      allowsnull=False,
                      help="Every X amount of steps of conjugate gradient algorithm it will be performed and steepest "
                           "descent algorithm.")

        form.addParam('EMAlgorithmSystem', params.EnumParam,
                      label='Algorithm for system energy minimization',
                      expertLevel=LEVEL_ADVANCED,
                      display=EnumParam.DISPLAY_HLIST,
                      choices=self.EM_PARAMS_OPTIONS,
                      default=1,
                      allowsnull=False,
                      help='Algorithm used to run the energy minimization for the hydrogens.\n*"steep"* will use the '
                           'steepest descent algorithm.\n*"cg"* will use a conjugate gradient algorithm which will also'
                           ' do steepest descent every X steeps.\nSteepest descent is more accurate and time consuming'
                           ' than conjugate gradient. It is recommended to perform steepest descent for the first 2 '
                           'minimizations and conjugate gradient for the last one. In case you only perform 1 '
                           'minimzation it is recommended to perform steepest descent which is more accurate.')

        form.addParam('EMNstcgsteepSystem', params.IntParam,
                      label='Nstcgsteep for CG algorithm in System minimization:',
                      expertLevel=LEVEL_ADVANCED,
                      condition='EMAlgorithmSystem == 1',
                      default=1000,
                      allowsnull=False,
                      help="Every X amount of steps of conjugate gradient algorithm it will be performed and steepest "
                           "descent algorithm.")

        form.addParam('ParamsFile', params.PathParam,
                      label='Parameters file',
                      condition="ScipionOrLocalParams == False",
                      allowsnull=False,
                      help="Input parameters file which will be used to run energy minimization")

        form.addParam('nsteps', params.IntParam,
                      label='Number of steps (nsteps)', default=50000,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Maximum number of steps to integrate or minimize, -1 is no maximum. "
                           "Multiplied by dt, it gives the ps of the simulation. e.g. 50000 nsteps * 0.002 dt = 100 ps."
                           " It depends on the system, but for NVT equilibration, 100 ps should be enough")

        form.addParam('UseGroFiles', params.PointerParam, label="Input Gro File",
                      pointerClass='GroFiles',
                      allowsNull=True,
                      help='Gro file which will be used to do energy minimization')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        if self.EMSteeps == 0:
            self._insertFunctionStep('getGromppParams')
        else:
            pass
        if self.EMSteeps == 1:
            self._insertFunctionStep('getGromppParams2')
        else:
            pass
        if self.EMSteeps == 2:
            self._insertFunctionStep('getGromppParams3')
        else:
            pass
        self._insertFunctionStep('createOutputStep')

    def getGromppParams(self):
        UseGroFile = os.path.abspath(self.UseGroFiles.get().getGro())
        UseTopolFile = os.path.abspath(self.UseGroFiles.get().getTopol())
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        nsteps = self.nsteps.get()
        working_dir = self.getWorkingDir()
        algorithm_hydrogen = self.EM_PARAMS_OPTIONS[self.EMAlgorithmHydrogen.get()]
        algorithm_water = self.EM_PARAMS_OPTIONS[self.EMAlgorithmWater.get()]
        algorithm_system = self.EM_PARAMS_OPTIONS[self.EMAlgorithmSystem.get()]
        if self.EMAlgorithmHydrogen.get() == 1:
            CG_steps_hydrogen = self.EMNstcgsteepHydrogen.get()
        else:
            CG_steps_hydrogen = 1000
        if self.EMAlgorithmWater.get() == 1:
            CG_steps_water = self.EMNstcgsteepWater.get()
        else:
            CG_steps_water = 1000
        if self.EMAlgorithmSystem.get() == 1:
            CG_steps_system = self.EMNstcgsteepSystem.get()
        else:
            CG_steps_system = 1000
        if self.ScipionOrLocalParams == True:
            minim_1_mdp = '%s/minim_1.mdp' % (working_dir)
            minim_1_mdp_file = "define = -DPOSRES_WATER -DPOSRES \n" \
                        "integrator = %s \n" \
                        "emtol = 1000.0 \n" \
                        "emstep = 0.01 \n" \
                        "nsteps = %d \n" \
                        "nstcgsteep = %d \n" \
                        "\n" \
                        "nstlist = 1 \n" \
                        "cutoff-scheme = Verlet \n" \
                        "ns_type = grid \n" \
                        "coulombtype = PME \n" \
                        "rcoulomb = 1.0 \n" \
                        "rvdw = 1.0 \n" \
                        "pbc = xyz" % (algorithm_hydrogen, nsteps, CG_steps_hydrogen)

            f = open(minim_1_mdp, "w")
            f.write(minim_1_mdp_file)
            f.close()
            mdp_file_1 = os.path.basename(minim_1_mdp)

            minim_2_mdp = '%s/minim_2.mdp' % (working_dir)
            minim_2_mdp_file = "define = -DPOSRES \n" \
                        "integrator = %s \n" \
                        "emtol = 1000.0 \n" \
                        "emstep = 0.01 \n" \
                        "nsteps = %d \n" \
                        "nstcgsteep = %d \n" \
                        "\n" \
                        "nstlist = 1 \n" \
                        "cutoff-scheme = Verlet \n" \
                        "ns_type = grid \n" \
                        "coulombtype = PME \n" \
                        "rcoulomb = 1.0 \n" \
                        "rvdw = 1.0 \n" \
                        "pbc = xyz" % (algorithm_water,nsteps, CG_steps_water)

            f = open(minim_2_mdp, "w")
            f.write(minim_2_mdp_file)
            f.close()
            mdp_file_2 = os.path.basename(minim_2_mdp)

            minim_3_mdp = '%s/minim_3.mdp' % (working_dir)
            minim_3_mdp_file = "define = -DFLEXIBLE \n" \
                        "integrator = %s \n" \
                        "emtol = 1000.0 \n" \
                        "emstep = 0.01 \n" \
                        "nsteps = %d \n" \
                        "nstcgsteep = %d \n" \
                        "\n" \
                        "nstlist = 1 \n" \
                        "cutoff-scheme = Verlet \n" \
                        "ns_type = grid \n" \
                        "coulombtype = PME \n" \
                        "rcoulomb = 1.0 \n" \
                        "rvdw = 1.0 \n" \
                        "pbc = xyz" % (algorithm_system, nsteps, CG_steps_system)

            f = open(minim_3_mdp, "w")
            f.write(minim_3_mdp_file)
            f.close()
            mdp_file_3 = os.path.basename(minim_3_mdp)
        else:
            UseMdpFile = self.ParamsFile.get()
            minim_mdp = UseMdpFile
            mdp_file_1 = minim_mdp
        if self.ScipionOrLocalParams == True:
            em_number = "_1"
        else:
            em_number = ""
        if self.useGpu.get():
            gpu = " -nb gpu"
        else:
            gpu = ""
        params_grompp_1 = 'grompp -f %s -c %s -r %s -p ' \
                        '%s -o em%s.tpr' % (mdp_file_1, UseGroFile, UseGroFile, UseTopolFile, em_number)
        params_grompp_2 = 'grompp -f %s -c em_1.gro -r em_1.gro -p ' \
                        '%s -o em_2.tpr' % (mdp_file_2, UseTopolFile)
        params_grompp_3 = 'grompp -f %s -c em_2.gro -r em_2.gro -p ' \
                        '%s -o em.tpr' % (mdp_file_3, UseTopolFile)
        params_genion_1 = ' mdrun -v -deffnm em%s%s' % (em_number, gpu)
        params_genion_2 = ' mdrun -v -deffnm em_2%s' % (gpu)
        params_genion_3 = ' mdrun -v -deffnm em%s' % (gpu)
        self.runJob(program, params_grompp_1, cwd=self._getPath())
        self.runJob(program, params_genion_1, cwd=self._getPath())

        if self.ScipionOrLocalParams == True:
            self.runJob(program, params_grompp_2, cwd=self._getPath())
            self.runJob(program, params_genion_2, cwd=self._getPath())
            self.runJob(program, params_grompp_3, cwd=self._getPath())
            self.runJob(program, params_genion_3, cwd=self._getPath())
        else:
            pass

    def getGromppParams2(self):
        UseGroFile = os.path.abspath(self.UseGroFiles.get().getGro())
        UseTopolFile = os.path.abspath(self.UseGroFiles.get().getTopol())
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        nsteps = self.nsteps.get()
        working_dir = self.getWorkingDir()
        algorithm_water = self.EM_PARAMS_OPTIONS[self.EMAlgorithmWater.get()]
        algorithm_system = self.EM_PARAMS_OPTIONS[self.EMAlgorithmSystem.get()]
        if self.EMAlgorithmWater.get() == 1:
            CG_steps_water = self.EMNstcgsteepWater.get()
        else:
            CG_steps_water = 1000
        if self.EMAlgorithmSystem.get() == 1:
            CG_steps_system = self.EMNstcgsteepSystem.get()
        else:
            CG_steps_system = 1000
        minim_2_mdp = '%s/minim_2.mdp' % (working_dir)
        minim_2_mdp_file = "define = -DPOSRES \n" \
                    "integrator = %s \n" \
                    "emtol = 1000.0 \n" \
                    "emstep = 0.01 \n" \
                    "nsteps = %d \n" \
                    "nstcgsteep = %d \n" \
                    "\n" \
                    "nstlist = 1 \n" \
                    "cutoff-scheme = Verlet \n" \
                    "ns_type = grid \n" \
                    "coulombtype = PME \n" \
                    "rcoulomb = 1.0 \n" \
                    "rvdw = 1.0 \n" \
                    "pbc = xyz" % (algorithm_water,nsteps, CG_steps_water)
        f = open(minim_2_mdp, "w")
        f.write(minim_2_mdp_file)
        f.close()
        mdp_file_2 = os.path.basename(minim_2_mdp)
        minim_3_mdp = '%s/minim_3.mdp' % (working_dir)
        minim_3_mdp_file = "define = -DFLEXIBLE \n" \
                    "integrator = %s \n" \
                    "emtol = 1000.0 \n" \
                    "emstep = 0.01 \n" \
                    "nsteps = %d \n" \
                    "nstcgsteep = %d \n" \
                    "\n" \
                    "nstlist = 1 \n" \
                    "cutoff-scheme = Verlet \n" \
                    "ns_type = grid \n" \
                    "coulombtype = PME \n" \
                    "rcoulomb = 1.0 \n" \
                    "rvdw = 1.0 \n" \
                    "pbc = xyz" % (algorithm_system, nsteps, CG_steps_system)
        f = open(minim_3_mdp, "w")
        f.write(minim_3_mdp_file)
        f.close()
        mdp_file_3 = os.path.basename(minim_3_mdp)
        if self.useGpu.get():
            gpu = " -nb gpu"
        else:
            gpu = ""
        params_grompp_2 = 'grompp -f %s -c %s -r %s -p ' \
                        '%s -o em_2.tpr' % (mdp_file_2, UseGroFile, UseGroFile, UseTopolFile)
        params_grompp_3 = 'grompp -f %s -c em_2.gro -r em_2.gro -p ' \
                        '%s -o em.tpr' % (mdp_file_3, UseTopolFile)
        params_genion_2 = ' mdrun -v -deffnm em_2%s' % (gpu)
        params_genion_3 = ' mdrun -v -deffnm em%s' % (gpu)
        self.runJob(program, params_grompp_2, cwd=self._getPath())
        self.runJob(program, params_genion_2, cwd=self._getPath())
        self.runJob(program, params_grompp_3, cwd=self._getPath())
        self.runJob(program, params_genion_3, cwd=self._getPath())

    def getGromppParams3(self):
        UseGroFile = os.path.abspath(self.UseGroFiles.get().getGro())
        UseTopolFile = os.path.abspath(self.UseGroFiles.get().getTopol())
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        nsteps = self.nsteps.get()
        working_dir = self.getWorkingDir()
        algorithm_system = self.EM_PARAMS_OPTIONS[self.EMAlgorithmSystem.get()]

        if self.EMAlgorithmSystem.get() == 1:
            CG_steps_system = self.EMNstcgsteepSystem.get()
        else:
            CG_steps_system = 1000
        minim_3_mdp = '%s/minim_3.mdp' % (working_dir)
        minim_3_mdp_file = "define = -DFLEXIBLE \n" \
                    "integrator = %s \n" \
                    "emtol = 1000.0 \n" \
                    "emstep = 0.01 \n" \
                    "nsteps = %d \n" \
                    "nstcgsteep = %d \n" \
                    "\n" \
                    "nstlist = 1 \n" \
                    "cutoff-scheme = Verlet \n" \
                    "ns_type = grid \n" \
                    "coulombtype = PME \n" \
                    "rcoulomb = 1.0 \n" \
                    "rvdw = 1.0 \n" \
                    "pbc = xyz" % (algorithm_system, nsteps, CG_steps_system)
        f = open(minim_3_mdp, "w")
        f.write(minim_3_mdp_file)
        f.close()
        mdp_file_3 = os.path.basename(minim_3_mdp)
        if self.useGpu.get():
            gpu = " -nb gpu"
        else:
            gpu = ""
        params_grompp_3 = 'grompp -f %s -c %s -r %s -p ' \
                        '%s -o em.tpr' % (mdp_file_3, UseGroFile, UseGroFile, UseTopolFile)
        params_genion_3 = ' mdrun -v -deffnm em%s' % (gpu)
        self.runJob(program, params_grompp_3, cwd=self._getPath())
        self.runJob(program, params_genion_3, cwd=self._getPath())

    def createOutputStep(self):

        em_edr_baseName = 'em.edr'
        em_gro_baseName = 'em.gro'
        em_log_baseName = 'em.log'
        em_tpr_baseName = 'em.tpr'
        em_trr_baseName = 'em.trr'
        FileName = 'EM_files'

        em_edr_localPath = relpath(abspath(self._getPath(em_edr_baseName)))
        em_gro_localPath = relpath(abspath(self._getPath(em_gro_baseName)))
        em_log_localPath = relpath(abspath(self._getPath(em_log_baseName)))
        em_tpr_localPath = relpath(abspath(self._getPath(em_tpr_baseName)))
        em_trr_localPath = relpath(abspath(self._getPath(em_trr_baseName)))

        em_object = grobj.EMgroFile(em_edr=em_edr_localPath,
                                       em_gro=em_gro_localPath,
                                       em_log=em_log_localPath,
                                       em_trr=em_trr_localPath,
                                       em_tpr=em_tpr_localPath)

        self._defineOutputs(outputEM=em_object)
        self._defineSourceRelation(self.UseGroFiles, em_object)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            if self.EMSteeps == 0:
                summary.append('This protocol has performed *3* energy minimization steps.\n'
                               'It has performed first an hydrogen minimization step, then a water molecule '
                               'minimization step and finally a system minimization step.')
            elif self.EMSteeps == 1:
                summary.append('This protocol has performed *2* energy minimization steps.\n'
                               'It has performed first an hydrogen and water molecule minimization step and finally a '
                               'system minimization step')
            elif self.EMSteeps == 2:
                summary.append('This protocol has performed *1* energy minimization step.\n'
                               'It has performed a system minimization step.')
            else:
                pass
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append('The methods used to perform the energy minimization protocol have been *"gmx grompp"* to '
                           'create the tpr files and *"gmx mdrun"* to run the minimization simulation.')

        return methods

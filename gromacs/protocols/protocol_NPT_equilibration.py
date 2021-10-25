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
This module will perform an equilibration NPT
"""
import os
from os.path import exists, basename, abspath, relpath, join
import pyworkflow.utils as pwutils
from pyworkflow.protocol import Protocol, params, Integer, EnumParam, StringParam, FloatParam
from pyworkflow.utils import Message, runJob, createLink
import pwem.objects as emobj
from gromacs.objects import *
import gromacs.objects as grobj
from pyworkflow.protocol.params import (LEVEL_ADVANCED, LEVEL_NORMAL, USE_GPU)
from pwem.protocols import EMProtocol

GROMACS_SPC = 0
GROMACS_SPCE = 1
GROMACS_TIP3P = 2


GROMACS_WATERFF_NAME = dict()
GROMACS_WATERFF_NAME[GROMACS_SPC] = 'steep'
GROMACS_WATERFF_NAME[GROMACS_SPCE] = 'mp'
GROMACS_WATERFF_NAME[GROMACS_TIP3P] = 'md-vv'


INTEGRATOR_LIST = [GROMACS_WATERFF_NAME[GROMACS_SPC], GROMACS_WATERFF_NAME[GROMACS_SPCE],
GROMACS_WATERFF_NAME[GROMACS_TIP3P]]

class GromacsNPTEquilibration(EMProtocol):
    """
    This protocol will perform NPT equilibration in orther to stabilize the density of the system.
    It has to be performed after NVT equilibration. It is also called isothermal-isobaric ensemble
    It is the previous step before Molecular Dynamics simulation.
    """
    _label = 'NPT equilibration'
    LINCS_ORDER_OPTIONS = [1, 2]
    PME_ORDER_OPTIONS = [4, 6, 8, 10]
    TCOUPL_OPTIONS = ['V-rescale', 'berendsen']
    GENVEL_OPTIONS = ['yes', 'no']
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
                      expertLevel=LEVEL_NORMAL,
                      condition="expertLevel == LEVEL_NORMAL",
                      default=True,
                      help="'Yes' = It uses the default params for NPT equilibration: \n"
                           "define = -DPOSRES ; position restrain the protein\n" \
                           "integrator = md ; leap-frog integrator\n" \
                           "nsteps = 50000 ; 2 * 50000 = 100 ps\n" \
                           "dt = 0.002 ; 2 fs\n" \
                           "nstxout = 500 ; save coordinates every 1.0 ps\n" \
                           "nstvout = 500 ; save velocities every 1.0 ps\n" \
                           "nstenergy = 500 ; save energies every 1.0 ps\n" \
                           "nstlog = 500 ; update log file every 1.0 ps\n" \
                           "continuation = yes ; restarting after NVT\n" \
                           "constraint_algorithm = lincs ; holonomic constraints\n" \
                           "constraints = h-bonds ; bonds involving H are constrained\n" \
                           "lincs_iter = 1 ; accuracy of LINCS\n" \
                           "lincs_order = 4 ; also related to accuracy\n" \
                           "cutoff-scheme = Verlet ; Buffered neighbor searching\n" \
                           "ns_type = grid ; search neighboring grid cells\n" \
                           "nstlist = 10 ; 20 fs, largerly irrelevant with Verlet\n" \
                           "rcoulomb = 1.0 ; short-range electrostatic cutoff (in nm)\n" \
                           "rvdw = 1.0 ; short-ragen van der Waals cutoff (in nm)\n" \
                           "DispCorr = EnerPres ; account for cut-off vdW scheme\n" \
                           "coulombtype = PME ; Particle Mesh Ewals for long-range electrostatics\n" \
                           "pme_order = 4 ; cubic interpolation\n" \
                           "fourierspacing = 0.16 ; grid spacing for FFT\n" \
                           "tcoupl = V-rescale ; modified Berendsen thermostat\n" \
                           "tc-grps = Protein Non-Protein ; two coupling groups - more accurate\n" \
                           "tau_t = 0.1    0.1 ; time constant, in ps\n" \
                           "ref_t = 300    300 ; reference temeprature, one for each group, in K\n" \
                           "pcoupl = Parrinello-Rahman ; Preassure coupling on in NPT\n"
                           "pcoupltype = isotropic ; uniform scaling of box vectors\n"
                           "tau_p = 2.0 ; time constant, in ps\n"
                           "ref_p = 1.0 ; reference pressure, in bar\n"
                           "compressibility = 4.5e-5 ; isothermal compressibility of water, bar⁻1\n"
                           "refcoord_scaling = com\n" \
                           "pbc = xyz ; 3-D periodic boundary conditions\n" \
                           "gen_vel = yes ; assign velocities from Maxwell distribution\n" \
                      "'No' == It uses the params from an input parameters file.")

        form.addParam('ParamsFile', params.PathParam,
                      label='Parameters file',
                      expertLevel=LEVEL_NORMAL,
                      condition="expertLevel == LEVEL_NORMAL and ScipionOrLocalParams == False",
                      allowsnull=False,
                      help="Input parameters file which will be used to run NPT equilibration")

        form.addParam('UseGroFiles', params.PointerParam, label="Input Set of Gromacs Files",
                      pointerClass='GroFiles',
                      allowsNull=True,
                      help='Input Topol File, Posre file and Gro file.')

        form.addParam('UseNVTGroFile', params.PointerParam, label="Input nvt.gro File",
                      pointerClass='NVTgroFile',
                      allowsNull=True,
                      important=True,
                      help='This nvt.gro file will be used to do NPT equilibration')

        form.addParam('rcoulomb', params.FloatParam,
                       label="Coulomb distance cut-off (rcoulomb)", default=1.0,
                       expertLevel=LEVEL_ADVANCED,
                       help="Distance for the Coulomb cut-off. It should be between 1.0 and 1.5")

        form.addParam('nsteps', params.IntParam,
                      label='Number of steps (nsteps)', default=50000,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Maximum number of steps to integrate or minimize, -1 is no maximum. "
                           "Multiplied by dt, it gives the ps of the simulation. e.g. 50000 nsteps * 0.002 dt = 100 ps."
                           " It depends on the system, but for NPT equilibration, 100 ps should be enough")

        form.addParam('dt', params.FloatParam,
                       label="Time step for integration (dt)", default=0.001,
                       expertLevel=LEVEL_ADVANCED,
                       help="Time step for integration. Multiplied by nsteps, it gives the ps of the simulation."
                            "e.g. 50000 nsteps * 0.002 dt = 50 ps. It depends on the system, but for NPT equilibration"
                            "100 ps should be enough")

        form.addParam('lincsIter', params.EnumParam,
                      label='Lincs iteration number (lincs_iter).',
                      choices=self.LINCS_ORDER_OPTIONS,
                      display=EnumParam.DISPLAY_HLIST,
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Number of iterations to correct for rotational lengthening in LINCS. For normal runs a"
                           " single step is sufficient, but for NVE runs where you want to conserve energy accurately "
                           "or for accurate energy minimization you might want to increase it to 2.")

        form.addParam('pmeOrder', params.EnumParam,
                      label='Interpolation order for PME (pme_order).',
                      choices=self.PME_ORDER_OPTIONS,
                      display=EnumParam.DISPLAY_HLIST,
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Interpolation order for PME (Particle-Mesh Ewald). 4 equals cubic interpolation. You might"
                           " try 6/8/10 when running in parallel and simultaneously decrease grid dimension")

        form.addParam('tcoupl', params.EnumParam,
                      label='Temperature coupling (tcoupl).',
                      choices=self.TCOUPL_OPTIONS,
                      display=EnumParam.DISPLAY_HLIST,
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="V-rescale: Temperature coupling using velocity rescaling with a stochastic term. This"
                           " thermostat is similar to Berendsen coupling, with the same scaling using tau-t, but the "
                           "stochastic term ensures that a proper canonical ensemble is generated. \n"
                           "Berendsen: Temperature coupling with a Berendsen thermostat to a bath with temperature "
                           "ref-t, with time constant tau-t. Several groups can be coupled seperately, these are "
                           "specified in the tc-grps field separated by spaces")

        form.addParam('refT', params.IntParam,
                      label='Reference temperature in Kelvin (ref_t)', default=300,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Reference temperature in Kelvin for the system.")

        form.addParam('genVel', params.EnumParam,
                      label='Velocity generation (gen_vel)',
                      choices=self.GENVEL_OPTIONS,
                      display=EnumParam.DISPLAY_HLIST,
                      default=1,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Velocity generation. 'Yes': Generate velocities in gmx grompp according to a Maxwell "
                           "distribution. 'No': Do not generate velocities. The velocities are set to zero when there "
                           "are no velocites in the input structure file.")

        form.addParam('genSeed', params.IntParam,
                      label='Generation seed (gen_seed)', default=-1,
                      expertLevel=LEVEL_ADVANCED,
                      allowsnull=False,
                      help="Seed used to initialize random generator for random velocities. When gen-seed is set to -1,"
                           " a pseudo random seed is used.")



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps

        self._insertFunctionStep('getGromppParams')
        self._insertFunctionStep('getGenionParams')
        self._insertFunctionStep('getGromppParams2')
        self._insertFunctionStep('getGenionParams2')
        self._insertFunctionStep('createOutputStep')

    def getGromppParams(self):
        UseTopolFile = os.path.abspath(self.UseGroFiles.get().getTopol())
        UsePosreFile = os.path.abspath(self.UseGroFiles.get().getPosre())
        UseNVTGroFile = os.path.abspath(self.UseNVTGroFile.get().getNVT_gro())
        UseCPTFile = os.path.abspath(self.UseNVTGroFile.get().getNVT_cpt())

        program = os.path.join("", '/usr/local/gromacs/bin/gmx')


        nsteps = self.nsteps.get()
        lincs_iter = self.LINCS_ORDER_OPTIONS[self.lincsIter.get()]
        rcoulomb = self.rcoulomb.get()
        pme_order = self.PME_ORDER_OPTIONS[self.pmeOrder.get()]
        if self.tcoupl.get() == 0:
            tcoupl = 'V-rescale'
        else:
            tcoupl = 'berendsen'
        ref_t = self.refT.get()
        if self.genVel.get() == 0:
            gen_vel = 'yes'
        else:
            gen_vel = 'no'
        gen_seed = self.genSeed.get()
        dt = self.dt.get()
        working_dir = self.getWorkingDir()
        if self.ScipionOrLocalParams == True:
            npt_mdp = '%s/npt.mdp' % (working_dir)
            npt_mdp_file ="title = OPLS Lysozyme NPT equilibration \n" \
                            "define = -DPOSRES \n" \
                            "integrator = md \n" \
                            "nsteps = %d \n" \
                            "dt = %f \n" \
                            "nstxout = 500 \n" \
                            "nstvout = 500 \n" \
                            "nstenergy = 500 \n" \
                            "nstlog = 500 \n" \
                            "continuation = yes \n" \
                            "constraint_algorithm = lincs \n" \
                            "constraints = h-bonds \n" \
                            "lincs_iter = %d \n" \
                            "lincs_order = 4 \n" \
                            "cutoff-scheme = Verlet \n" \
                            "ns_type = grid \n" \
                            "nstlist = 10 \n" \
                            "rcoulomb = %d \n" \
                            "rvdw = 1.0 \n" \
                            "DispCorr = EnerPres \n" \
                            "coulombtype = PME \n" \
                            "pme_order = %d \n" \
                            "fourierspacing = 0.16 \n" \
                            "tcoupl = %s \n" \
                            "tc-grps = Protein Non-Protein \n" \
                            "tau_t = 0.1    0.1 \n" \
                            "ref_t = %d    %d \n" \
                            "pcoupl = Parrinello-Rahman \n" \
                            "pcoupltype = isotropic \n" \
                            "tau_p = 2 \n" \
                            "ref_p = 1 \n" \
                            "compressibility = 4.5e-5 \n" \
                            "refcoord_scaling = com \n" \
                            "pbc = xyz \n" \
                            "gen_vel = %s \n" % (nsteps, dt, lincs_iter, rcoulomb, pme_order,
                                               tcoupl, ref_t, ref_t, gen_vel)

            f = open(npt_mdp, "w")
            f.write(npt_mdp_file)
            f.close()
            mdp_file = os.path.basename(npt_mdp)
        else:
            UseMdpFile = self.ParamsFile.get()
            npt_mdp = UseMdpFile
            mdp_file = npt_mdp

        params_grompp = 'grompp -f %s -c %s -r %s -t %s -p ' \
                        '%s -o npt2.tpr' % (mdp_file, UseNVTGroFile, UseNVTGroFile, UseCPTFile, UseTopolFile)
        self.runJob(program, params_grompp, cwd=self._getPath())

    def getGenionParams(self):
        if self.useGpu.get():
            gpu = " -nb gpu"
        else:
            gpu = ""
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        params_genion = ' mdrun -v -deffnm npt2%s' % (gpu)

        self.runJob(program, params_genion, cwd=self._getPath())

    def getGromppParams2(self):
        UseTopolFile = os.path.abspath(self.UseGroFiles.get().getTopol())
        UsePosreFile = os.path.abspath(self.UseGroFiles.get().getPosre())
        UseNVTGroFile = os.path.abspath(self.UseNVTGroFile.get().getNVT_gro())
        UseCPTFile = os.path.abspath(self.UseNVTGroFile.get().getNVT_cpt())

        program = os.path.join("", '/usr/local/gromacs/bin/gmx')

        nsteps = self.nsteps.get()
        lincs_iter = self.LINCS_ORDER_OPTIONS[self.lincsIter.get()]
        rcoulomb = self.rcoulomb.get()
        pme_order = self.PME_ORDER_OPTIONS[self.pmeOrder.get()]
        if self.tcoupl.get() == 0:
            tcoupl = 'V-rescale'
        else:
            tcoupl = 'berendsen'
        ref_t = self.refT.get()
        if self.genVel.get() == 0:
            gen_vel = 'yes'
        else:
            gen_vel = 'no'
        gen_seed = self.genSeed.get()
        dt = self.dt.get()
        working_dir = self.getWorkingDir()
        if self.ScipionOrLocalParams == True:
            npt_mdp_2 = '%s/npt_2.mdp' % (working_dir)
            npt_mdp_file_2 ="title = OPLS Lysozyme NPT equilibration \n" \
                            "define = -DPOSRES_LOW \n" \
                            "integrator = md \n" \
                            "nsteps = %d \n" \
                            "dt = %f \n" \
                            "nstxout = 500 \n" \
                            "nstvout = 500 \n" \
                            "nstenergy = 500 \n" \
                            "nstlog = 500 \n" \
                            "continuation = yes \n" \
                            "constraint_algorithm = lincs \n" \
                            "constraints = h-bonds \n" \
                            "lincs_iter = %d \n" \
                            "lincs_order = 4 \n" \
                            "cutoff-scheme = Verlet \n" \
                            "ns_type = grid \n" \
                            "nstlist = 10 \n" \
                            "rcoulomb = %d \n" \
                            "rvdw = 1.0 \n" \
                            "DispCorr = EnerPres \n" \
                            "coulombtype = PME \n" \
                            "pme_order = %d \n" \
                            "fourierspacing = 0.16 \n" \
                            "tcoupl = %s \n" \
                            "tc-grps = Protein Non-Protein \n" \
                            "tau_t = 0.1    0.1 \n" \
                            "ref_t = %d    %d \n" \
                            "pcoupl = Parrinello-Rahman \n" \
                            "pcoupltype = isotropic \n" \
                            "tau_p = 2 \n" \
                            "ref_p = 1 \n" \
                            "compressibility = 4.5e-5 \n" \
                            "refcoord_scaling = com \n" \
                            "pbc = xyz \n" \
                            "gen_vel = %s \n" % (nsteps, dt, lincs_iter, rcoulomb, pme_order,
                                               tcoupl, ref_t, ref_t, gen_vel)
            f = open(npt_mdp_2, "w")
            f.write(npt_mdp_file_2)
            f.close()
            mdp_file_2 = os.path.basename(npt_mdp_2)
        else:
            UseMdpFile = self.ParamsFile.get()
            npt_mdp_2 = UseMdpFile
            mdp_file_2 = npt_mdp_2

        params_grompp = 'grompp -f %s -c npt2.gro -r npt2.gro -t npt2.cpt -p ' \
                        '%s -o npt.tpr' % (mdp_file_2, UseTopolFile)
        self.runJob(program, params_grompp, cwd=self._getPath())

    def getGenionParams2(self):
        if self.useGpu.get():
            gpu = " -nb gpu"
        else:
            gpu = ""
        program = os.path.join("", '/usr/local/gromacs/bin/gmx')
        params_genion = ' mdrun -v -deffnm npt%s' % (gpu)

        self.runJob(program, params_genion, cwd=self._getPath())

    def createOutputStep(self):
        npt_cpt_baseName = 'npt.cpt'
        npt_edr_baseName = 'npt.edr'
        npt_gro_baseName = 'npt.gro'
        npt_log_baseName = 'npt.log'
        npt_tpr_baseName = 'npt.tpr'
        npt_trr_baseName = 'npt.trr'
        npt_cpt_localPath = relpath(abspath(self._getPath(npt_cpt_baseName)))
        npt_edr_localPath = relpath(abspath(self._getPath(npt_edr_baseName)))
        npt_gro_localPath = relpath(abspath(self._getPath(npt_gro_baseName)))
        npt_log_localPath = relpath(abspath(self._getPath(npt_log_baseName)))
        npt_tpr_localPath = relpath(abspath(self._getPath(npt_tpr_baseName)))
        npt_trr_localPath = relpath(abspath(self._getPath(npt_trr_baseName)))
        npt_object = grobj.NPTgroFile(npt_edr=npt_edr_localPath,
                                       npt_gro=npt_gro_localPath,
                                       npt_log=npt_log_localPath,
                                       npt_trr=npt_trr_localPath,
                                       npt_tpr=npt_tpr_localPath,
                                       npt_cpt=npt_cpt_localPath)
        self._defineOutputs(outputNPT=npt_object)
        self._defineSourceRelation(self.UseGroFiles, npt_object)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        nsteps = self.nsteps.get()
        dt = self.dt.get()
        time = nsteps*dt
        if self.isFinished():
            summary.append('This protocol has performed a NPT equilibration satisfactorily'
                           '\nThe equilibration NPT has been performed during'
                           ' *%d ps*' % (time))
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append('The methods used to perform the NPT equilibration protocol have been *"gmx grompp"* to '
                           'create the tpr files and *"gmx mdrun"* to run the NPT equilibration.')

        return methods

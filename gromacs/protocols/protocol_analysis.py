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
This module will perform an analysis of the Molecular Dynamics
"""
import os
from os.path import exists, basename, abspath, relpath, join
import pyworkflow.utils as pwutils
from pyworkflow.protocol import Protocol, params, Integer, EnumParam, StringParam, FloatParam
from pyworkflow.utils import Message, runJob, createLink
import pwem.objects as emobj
from gromacs.objects import *
import gromacs.objects as grobj
from pyworkflow.protocol.params import (LEVEL_ADVANCED, USE_GPU)
from pwem.protocols import EMProtocol

class GromacsAnalysis(EMProtocol):
    """
    This protocol will perform a Molecular Dynamics analysis.
    It takes the files from the "Molecular Dynamics" protocol.
    It is recommended to perform the trjconv analysis in order to fix the trajectory.
    """
    _label = 'MD Analysis'
    LINCS_ORDER_OPTIONS = [1, 2]
    PME_ORDER_OPTIONS = [4, 6, 8, 10]
    TCOUPL_OPTIONS = ['V-rescale', 'berendsen']
    GENVEL_OPTIONS = ['yes', 'no']
    ANALYSIS_OPTIONS = ['trjconv', 'rms', 'rmsf', 'gyrate', 'sasa', 'hbond']
    XTC_SOURCE_OPTIONS = ['MD', 'Local file']
    # -------------------------- DEFINE constants ----------------------------


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """

        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('XtcFileSource', params.EnumParam,
                      label='Choose files source ',
                      choices=self.XTC_SOURCE_OPTIONS,
                      display=EnumParam.DISPLAY_HLIST,
                      default=0,
                      allowsnull=False,
                      help='Choose source for the xtc file. trjconv takes it from a previous trjconv analysis run.\n'
                           ' raw MD takes it from the molecular dynamics without preprocessing with trjconv.\n'
                           ' Local file takes a file from a local directory.')

        form.addParam('NdxAnalysis', params.BooleanParam,
                      label="Perform analysis on NDX selection",
                      default=False,
                      important=True,
                      help='Perform all analysis on an NDX selection of your own.')

        form.addParam('NdxFileLocal', params.PathParam,
                      label='Input ndx file from local folder',
                      allowsnull=False,
                      condition='NdxAnalysis == True',
                      help='Choose a local ndx file to perform the analysis.')

        form.addParam('NdxGroup', params.IntParam,
                      label="Select ndx group to analyse",
                      default=18,
                      condition='NdxAnalysis == True',
                      allowsnull=False,
                      helo="Group taken from ndx file to analyse.")

        form.addParam('TrjconvAnalysis', params.BooleanParam,
                      label="trjconv analysis",
                      default=True,
                      help='Perform trjconv analysis to fix the trajectory file. Sometimes, the atoms of the protein go'
                           'out the box and it has to be fixed. This protocol also takes the first frame as a'
                           'reference to put on it all the subsequent frames')

        form.addParam('RmsAnalysis', params.BooleanParam,
                      label="rms analysis",
                      default=True,
                      help='Perform RMS analysis on the protein or selection')

        form.addParam('RmsfAnalysis', params.BooleanParam,
                      label="rmsf analysis",
                      default=True,
                      help='Perform RMSF analysis on the protein or selection')

        form.addParam('SasaAnalysis', params.BooleanParam,
                      label="sasa analysis",
                      default=True,
                      help='Perform trjconv analysis to fix trajectory file.')

        form.addParam('HbondAnalysis', params.BooleanParam,
                      label="hbond analysis",
                      default=True,
                      help='Perform hbond analysis to visualyze how many atoms a residue has along the simulation.')

        form.addParam('ResidueHbond', params.IntParam,
                      label="Input residue to perform Hbond analysis",
                      allowsnull=False,
                      condition='HbondAnalysis == True',
                      help='Choose a residue number to perform Hbond analysis')

        form.addParam('GroFileLocal', params.PathParam,
                      label='Input gro file from local folder',
                      allowsnull=False,
                      condition='XtcFileSource == 1 and HbondAnalysis == True',
                      help='Choose a local tpr file to perform rms analysis.')

        form.addParam('ClusterAnalysis', params.BooleanParam,
                      label="cluster analysis",
                      default=True,
                      help='Perform cluster analysis to take the most representative frames and their structures')

        form.addParam('ClusterCutoff', params.FloatParam,
                      label="Enter cutoff for cluster analysis",
                      default=0.1,
                      condition='ClusterAnalysis == True',
                      allowsnull=False,
                      helo="Cutoff used for cluster analysis.")

        form.addParam('DensityAnalysis', params.BooleanParam,
                      label="Density analysis",
                      default=True,
                      help='Perform densisty analysis on the protein or selection.')

        form.addParam('GetReferenceStructure', params.BooleanParam,
                      label="Get Reference Structure from the first MD frame in pdb format",
                      default=True,
                      help='Perform trjcon command to take the reference structure of the simulation. It takes the'
                           'first frame.')

        form.addParam('ReduceFrames', params.BooleanParam,
                      label="Reduce the number of frames",
                      default=True,
                      help='Perform trjconv commad to reduce the number of frames by a certain number')

        form.addParam('dtFrames', params.IntParam,
                      label="Extract frames every dt =",
                      default=10,
                      condition="ReduceFrames == True",
                      help="Extract a portion of frames. e.g. 10 == Extract 1/10 frames")

        form.addParam('XtcFileLocal', params.PathParam,
                      label='Input xtc file from local folder',
                      allowsnull=False,
                      condition='XtcFileSource == 1',
                      help='Choose a local xtc file to perform the analysis.')

        form.addParam('TprFileLocal', params.PathParam,
                      label='Input tpr file from local folder',
                      allowsnull=False,
                      condition='XtcFileSource == 1',
                      help='Choose a local tpr file to perform the analysis.')

        form.addParam('UseMDGroFile', params.PointerParam, label="Input MD files from MD run",
                      pointerClass='MDgroFile',
                      allowsNull=False,
                      condition='XtcFileSource == 0',
                      help='This molecular dynamics set of files will be used to perform the analysis.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        if self.TrjconvAnalysis.get() == True:
            self._insertFunctionStep('getTrjconvParams')
        if self.RmsAnalysis.get() == True:
            self._insertFunctionStep('getRmsParams')
        if self.RmsfAnalysis.get() == True:
            self._insertFunctionStep('getRmsfParams')
        if self.SasaAnalysis.get() == True:
            self._insertFunctionStep('getSasaParams')
        if self.HbondAnalysis.get() == True:
            self._insertFunctionStep('getHbondParams')
        if self.ClusterAnalysis.get() == True:
            self._insertFunctionStep('getClusterParams')
        if self.DensityAnalysis.get() == True:
            self._insertFunctionStep('getDensityParams')
        if self.GetReferenceStructure.get() == True:
            self._insertFunctionStep('getReferenceStructureParams')
        if self.ReduceFrames.get() == True:
            self._insertFunctionStep('getReduceFramesParams')
        else:
            print("nothing")
        self._insertFunctionStep('xmgrace')
        self._insertFunctionStep('createOutputStep')

    def getTrjconvParams(self):
        if self.XtcFileSource == 0:
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        else:
            UseMDtprFile = os.path.abspath(self.TprFileLocal.get())
            UseMDxtcFile = os.path.abspath(self.XtcFileLocal.get())
            new_file = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
        program = 'printf "4\n0" | /usr/local/gromacs/bin/gmx'
        params_trjconv = 'trjconv -s %s -f %s -o %s_pre_noPBC.xtc ' \
                         '-pbc mol -center' % (UseMDtprFile, UseMDxtcFile,
                                               new_file)
        program_2 = 'printf "4\n0" | /usr/local/gromacs/bin/gmx'
        params_trjconv_2 = 'trjconv -s %s -f %s_pre_noPBC.xtc -fit progressive -o %s_noPBC.xtc' % (UseMDtprFile,
                                                                                                   new_file, new_file)
        self.runJob(program, params_trjconv, cwd=self._getPath())
        self.runJob(program_2, params_trjconv_2, cwd=self._getPath())

    def getRmsParams(self):
        if self.NdxAnalysis.get() == True:
            ndx_file = ' -n %s' % (self.NdxFileLocal.get())
            ndx_group = self.NdxGroup.get()
        else:
            ndx_file = ''
            ndx_group=4
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "%d\n%d" | /usr/local/gromacs/bin/gmx' % (ndx_group, ndx_group)
        params_rms = 'rms -s %s -f %s -o %s_rmsd.xvg ' \
                        '-tu ns%s' % (UseMDtprFile, UseMDxtcFile, new_file, ndx_file)

        self.runJob(program, params_rms, cwd=self._getPath())

    def getRmsfParams(self):
        if self.NdxAnalysis.get() == True:
            ndx_file = ' -n %s' % (self.NdxFileLocal.get())
            ndx_group = self.NdxGroup.get()
        else:
            ndx_file = ''
            ndx_group=1
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "%d" | /usr/local/gromacs/bin/gmx' % (ndx_group)
        params_rmsf = 'rmsf -s %s -f %s -o %s_rmsf.xvg ' \
                        '-res%s' % (UseMDtprFile, UseMDxtcFile, new_file, ndx_file)

        self.runJob(program, params_rmsf, cwd=self._getPath())

    def getSasaParams(self):
        if self.NdxAnalysis.get() == True:
            ndx_file = ' -n %s' % (self.NdxFileLocal.get())
            ndx_group = self.NdxGroup.get()
        else:
            ndx_file = ''
            ndx_group=1
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "%d" | /usr/local/gromacs/bin/gmx' % (ndx_group)
        params_sasa = 'sasa -s %s -f %s -o %s_sasa2.xvg -oa %s_atomic2-sasa.xvg -or %s_residue2-sasa.xvg%s' % (
            UseMDtprFile, UseMDxtcFile, new_file, new_file, new_file, ndx_file)

        self.runJob(program, params_sasa, cwd=self._getPath())

    def getHbondParams(self):
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
            UseMDgroFile = os.path.abspath(self.UseMDGroFile.get().getMD_gro())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
            UseMDgroFile = os.path.abspath(self.UseMDGroFile.get().getMD_gro())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
            UseMDgroFile = self.GroFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
            UseMDgroFile = self.GroFileLocal.get()
        residue = self.ResidueHbond.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program_ndx = 'printf "r%d\n1 &! r%d\nq\n" | /usr/local/gromacs/bin/gmx' % (residue, residue)
        params_ndx = 'make_ndx -f %s -o %s_hbond.ndx' % (UseMDgroFile, new_file)
        program = 'printf "Protein_&_!r_%s\nr_%s" | /usr/local/gromacs/bin/gmx' % (residue, residue)
        params_hbond_1 = 'hbond -s %s -f %s -num %s_hydrogen2-bond-intra-protein.xvg -n %s_hbond.ndx' % (UseMDtprFile, UseMDxtcFile,
                                                                                      new_file, new_file)
        self.runJob(program_ndx, params_ndx, cwd=self._getPath())
        self.runJob(program, params_hbond_1, cwd=self._getPath())

    def getClusterParams(self):
        if self.NdxAnalysis.get() == True:
            ndx_file = ' -n %s' % (self.NdxFileLocal.get())
            ndx_group = self.NdxGroup.get()
        else:
            ndx_file = ''
            ndx_group=4
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        cutoff = self.ClusterCutoff.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "%d\n1" | /usr/local/gromacs/bin/gmx' % (ndx_group)
        params_cluster = 'cluster -s %s -f %s -method gromos -cl %s_average.pdb -g %s_cluster.log ' \
                      '-cutoff %f%s -av' % (UseMDtprFile, UseMDxtcFile, new_file, new_file, cutoff, ndx_file)

        self.runJob(program, params_cluster, cwd=self._getPath())

    def getDensityParams(self):
        if self.NdxAnalysis.get() == True:
            ndx_file = ' -n %s' % (self.NdxFileLocal.get())
            ndx_group = self.NdxGroup.get()
        else:
            ndx_file = ''
            ndx_group=1
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "%d" | /usr/local/gromacs/bin/gmx' % (ndx_group)
        params_density = 'density -s %s -f %s -o %s_density.xvg%s ' % (UseMDtprFile, UseMDxtcFile, new_file, ndx_file)
        self.runJob(program, params_density, cwd=self._getPath())

    def getReduceFramesParams(self):
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        dt_frames = self.dtFrames.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "Protein" | /usr/local/gromacs/bin/gmx'
        params_reduce = 'trjconv -s %s -f %s  -o %s_reduced_traj.pdb -dt %d' % (UseMDtprFile, UseMDxtcFile,
                                                                                      new_file, dt_frames)
        self.runJob(program, params_reduce, cwd=self._getPath())

    def getReferenceStructureParams(self):
        if self.TrjconvAnalysis.get() == False and self.XtcFileSource == 0:
            UseMDxtcFile = os.path.abspath(self.UseMDGroFile.get().getMD_xtc())
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == True and self.XtcFileSource == 0:
            new_file_xtc = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = os.path.abspath(self.UseMDGroFile.get().getMD_tpr())
        elif self.TrjconvAnalysis.get() == False and self.XtcFileSource == 1:
            UseMDxtcFile = self.XtcFileLocal.get()
            UseMDtprFile = self.TprFileLocal.get()
        else:
            new_file_xtc = basename(os.path.abspath(self.XtcFileLocal.get()).split(".")[0])
            UseMDxtcFile = '%s_noPBC.xtc' % (new_file_xtc)
            UseMDtprFile = self.TprFileLocal.get()
        dt_frames = self.dtFrames.get()
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        program = 'printf "Protein" | /usr/local/gromacs/bin/gmx'
        params_reference = 'trjconv -s %s -f %s  -o %s_reference_structure.pdb -dump -0' % (UseMDtprFile, UseMDxtcFile,
                                                                                      new_file)
        self.runJob(program, params_reference, cwd=self._getPath())

    def xmgrace(self):
        if self.RmsAnalysis == True:
            if self.XtcFileSource == 1:
                new_file = basename((self.XtcFileLocal.get()).split(".")[0])
            else:
                new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            program = os.path.join("", 'grace')
            params_grace = '%s_rmsd.xvg -hdevice PNG -hardcopy -printfile %s_rmsd.png' % (new_file, new_file)
            self.runJob(program, params_grace, cwd=self._getPath())
        else:
            pass
        if self.RmsfAnalysis == True:
            if self.XtcFileSource == 1:
                new_file = basename((self.XtcFileLocal.get()).split(".")[0])
            else:
                new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            program = os.path.join("", 'grace')
            params_grace = '%s_rmsf.xvg -hdevice PNG -hardcopy -printfile %s_rmsf.png' % (new_file, new_file)
            self.runJob(program, params_grace, cwd=self._getPath())
        else:
            pass
        if self.SasaAnalysis == True:
            if self.XtcFileSource == 1:
                new_file = basename((self.XtcFileLocal.get()).split(".")[0])
            else:
                new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            program = os.path.join("", 'grace')
            params_grace_1 = '%s_sasa2.xvg -hdevice PNG -hardcopy -printfile %s_sasa2.png' % (new_file, new_file)
            params_grace_2 = '%s_atomic2-sasa.xvg -hdevice PNG -hardcopy -printfile ' \
                             '%s_atomic2-sasa.png' % (new_file, new_file)
            params_grace_3 = '%s_residue2-sasa.xvg -hdevice PNG -hardcopy -printfile ' \
                             '%s_residue2-sasa.png' % (new_file, new_file)
            self.runJob(program, params_grace_1, cwd=self._getPath())
            self.runJob(program, params_grace_2, cwd=self._getPath())
            self.runJob(program, params_grace_3, cwd=self._getPath())
        else:
            pass
        if self.HbondAnalysis == True:
            if self.XtcFileSource == 1:
                new_file = basename((self.XtcFileLocal.get()).split(".")[0])
            else:
                new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            program = os.path.join("", 'grace')
            params_grace_1 = '%s_hydrogen2-bond-intra-protein.xvg -hdevice PNG -hardcopy ' \
                             '-printfile %s_hydrogen2-bond-intra-protein.png' % (new_file, new_file)
            self.runJob(program, params_grace_1, cwd=self._getPath())
        else:
            pass
        if self.ClusterAnalysis == True:
            if self.XtcFileSource == 1:
                new_file = basename((self.XtcFileLocal.get()).split(".")[0])
            else:
                new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            program = os.path.join("", 'grace')
            params_grace_1 = 'rmsd-dist.xvg -hdevice PNG -hardcopy ' \
                             '-printfile rmsd-dist.png'
            self.runJob(program, params_grace_1, cwd=self._getPath())
        else:
            pass
        if self.DensityAnalysis == True:
            if self.XtcFileSource == 1:
                new_file = basename((self.XtcFileLocal.get()).split(".")[0])
            else:
                new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
            program = os.path.join("", 'grace')
            params_grace_1 = '%s_density.xvg -hdevice PNG -hardcopy -printfile %s_density.png' % (new_file, new_file)
            self.runJob(program, params_grace_1, cwd=self._getPath())
        else:
            pass

    def createOutputStep(self):
        if self.XtcFileSource == 1:
            new_file = basename((self.XtcFileLocal.get()).split(".")[0])
        else:
            new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
        if self.TrjconvAnalysis == True:
            md_0_1_xtc_noPDB_baseName = '%s_noPDB.xtc' % (new_file)
            md_0_1_xtc_noPDB_localPath = relpath(abspath(self._getPath(md_0_1_xtc_noPDB_baseName)))
            analysis_object_trjconv = grobj.TrjconvFile(trjconv=md_0_1_xtc_noPDB_localPath)
            self._defineOutputs(outputAnalysis_trjconv=analysis_object_trjconv)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_trjconv)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_trjconv)
        else:
            pass
        if self.RmsAnalysis == True:
            rmsd_xvg_baseName = '%s_rmsd.xvg' % (new_file)
            rmsd_xvg_localPath = relpath(abspath(self._getPath(rmsd_xvg_baseName)))
            analysis_object_rms = grobj.XvgFile(xvg=rmsd_xvg_localPath)
            self._defineOutputs(outputAnalysis_rms=analysis_object_rms)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_rms)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_rms)
        else:
            pass
        if self.RmsfAnalysis == True:
            rmsf_xvg_baseName = '%s_rmsf.xvg' % (new_file)
            rmsf_xvg_localPath = relpath(abspath(self._getPath(rmsf_xvg_baseName)))
            analysis_object_rmsf = grobj.XvgFile(xvg=rmsf_xvg_localPath)
            self._defineOutputs(outputAnalysis_rmsf=analysis_object_rmsf)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_rmsf)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_rmsf)
        else:
            pass
        if self.SasaAnalysis == True:
            sasa2_xvg_baseName = '%s_sasa2.xvg' % (new_file)
            sasa2_xvg_localPath = relpath(abspath(self._getPath(sasa2_xvg_baseName)))
            analysis_object_sasa_1 = grobj.XvgFile(xvg=sasa2_xvg_localPath)
            atomic2_sasa_xvg_baseName = '%s_atomic2-sasa.xvg' % (new_file)
            atomic2_sasa_xvg_localPath = relpath(abspath(self._getPath(atomic2_sasa_xvg_baseName)))
            analysis_object_sasa_2 = grobj.XvgFile(xvg=atomic2_sasa_xvg_localPath)
            residue2_sasa_xvg_baseName = '%s_residue2-sasa.xvg' % (new_file)
            residue2_sasa_xvg_localPath = relpath(abspath(self._getPath(residue2_sasa_xvg_baseName)))
            analysis_object_sasa_3 = grobj.XvgFile(xvg=residue2_sasa_xvg_localPath)
            self._defineOutputs(outputAnalysis_sasa1=analysis_object_sasa_1)
            self._defineOutputs(outputAnalysis_sasa2=analysis_object_sasa_2)
            self._defineOutputs(outputAnalysis_sasa3=analysis_object_sasa_3)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_sasa_1)
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_sasa_2)
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_sasa_3)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_sasa_1)
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_sasa_2)
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_sasa_3)
        else:
            pass
        if self.HbondAnalysis == True:
            hbond1_xvg_baseName = '%s_hydrogen2-bond-intra-protein.xvg' % (new_file)
            hbond1_xvg_localPath = relpath(abspath(self._getPath(hbond1_xvg_baseName)))
            analysis_object_hbond_1 = grobj.XvgFile(xvg=hbond1_xvg_localPath)
            self._defineOutputs(outputAnalysis_hbond1=analysis_object_hbond_1)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_hbond_1)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_hbond_1)
        else:
            pass
        if self.ClusterAnalysis == True:
            clusterpdb_xvg_baseName = '%s_average.pdb' % (new_file)
            clusterpdb_xvg_localPath = relpath(abspath(self._getPath(clusterpdb_xvg_baseName)))
            analysis_object_clusterpdb = grobj.XvgFile(xvg=clusterpdb_xvg_localPath)
            self._defineOutputs(outputAnalysis_clusterpdb=analysis_object_clusterpdb)
            clusterlog_xvg_baseName = '%s_cluster.log' % (new_file)
            clusterlog_xvg_localPath = relpath(abspath(self._getPath(clusterlog_xvg_baseName)))
            analysis_object_clusterlog = grobj.XvgFile(xvg=clusterlog_xvg_localPath)
            self._defineOutputs(outputAnalysis_clusterlog=analysis_object_clusterlog)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_clusterpdb)
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_clusterlog)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_clusterpdb)
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_clusterlog)
        else:
            pass

        if self.DensityAnalysis == True:
            density_xvg_baseName = '%s_density_protein.xvg' % (new_file)
            density_xvg_localPath = relpath(abspath(self._getPath(density_xvg_baseName)))
            analysis_object_density = grobj.XvgFile(xvg=density_xvg_localPath)
            self._defineOutputs(outputAnalysis_density=analysis_object_density)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_density)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_density)
        else:
            pass
        if self.GetReferenceStructure == True:
            pdb_reference_baseName = '%s_reference_structure.pdb' % (new_file)
            pdb_reference_localPath = relpath(abspath(self._getPath(pdb_reference_baseName)))
            analysis_object_pdb_reference = emobj.PdbFile(pdb=pdb_reference_localPath)
            self._defineOutputs(outputAnalysis_reference=analysis_object_pdb_reference)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_pdb_reference)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_pdb_reference)
        else:
            pass
        if self.ReduceFrames == True:
            pdb_reduced_baseName = '%s_reduced_traj.pdb' % (new_file)
            pdb_reduced_localPath = relpath(abspath(self._getPath(pdb_reduced_baseName)))
            analysis_object_pdb_reduced = emobj.PdbFile(pdb=pdb_reduced_localPath)
            self._defineOutputs(outputAnalysis_reduced=analysis_object_pdb_reduced)
            if self.XtcFileSource == 0:
                self._defineSourceRelation(self.UseMDGroFile, analysis_object_pdb_reduced)
            else:
                self._defineSourceRelation(self.XtcFileLocal, analysis_object_pdb_reduced)
        else:
            pass


    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []

        if self.isFinished():
            if self.TrjconvAnalysis:
                summary.append('This protocol has performed a Trjconv command, fixing the structure in case it has '
                               'atom position issues satisfactorily')
            if self.RmsAnalysis:
                summary.append('This protocol has performed RMS analysis')
            if self.RmsfAnalysis:
                summary.append('This protocol has performed RMSF analysis')
            if self.SasaAnalysis:
                summary.append('This protocol has performed Solvent accessible surface area analysis')
            if self.HbondAnalysis:
                summary.append('This protocol has performed Hbond analysis '
                               'on residue *%d*' % (self.ResidueHbond.get()))
            if self.ClusterAnalysis:
                summary.append('This protocol has performed Cluster analysis '
                               'with cutoff *%d*' % (self.ClusterCutoff.get()))
            if self.DensityAnalysis:
                summary.append('This protocol has performed Density analysis')
            if self.GetReferenceStructure:
                summary.append('This protocol has got the reference structure of the simulation')
            if self.ReduceFrames:
                summary.append('This protocol has reduced the frames by *%d*' % (self.dtFrames.get()))
        return summary

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append('This protocol has performed a complete molecular dynamics analysis satisfactorily')

        return methods

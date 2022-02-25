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

import os, glob, subprocess
import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params
from pwchem.viewers import VmdViewPopen
from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.viewers import ChimeraViewer

from pwchem.viewers import PyMolViewer, PyMolView
from pwchem.utils import natural_sort

from gromacs import Plugin
from ..objects import GromacsSystem
from ..protocols import GromacsMDSimulation
from ..constants import *

program = Plugin.getGromacsBin()

class GromacsSystemViewer(pwviewer.Viewer):
  _label = 'Viewer Gromacs system'
  _environments = [pwviewer.DESKTOP_TKINTER]
  _targets = [GromacsSystem]

  def _visualize(self, obj, **kwargs):
    groFile = os.path.abspath(obj.getSystemFile())

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(groFile, cwd=os.path.dirname(groFile))

class GromacsSimulationViewer(pwviewer.ProtocolViewer):
    """ Visualize the output of Desmond simulation """
    _label = 'Viewer Gromacs Simulation'
    _targets = [GromacsMDSimulation]
    _analysis = ['RMSD', 'RMSF', 'Gyration', 'SASA', 'HBond', 'Clustering']
    _ndxGroups = ['System', 'Protein', 'Protein-H', 'C-alpha', 'Backbone', 'MainChain',
                  'MainChain+Cb', 'MainChain+H', 'SideChain']

    def __init__(self, **args):
      super().__init__(**args)

    def _defineParams(self, form):
      form.addSection(label='Visualization of Gromacs Simulation')
      group = form.addGroup('Open Gromacs system')
      group.addParam('displayPymol', params.LabelParam,
                     label='Open system in PyMol: ',
                     help='Display System in Pymol GUI.'
                     )
      group = form.addGroup('Open MD simulation')
      group.addParam('chooseStage', params.EnumParam,
                     choices=self._getStagesWTrj(), default=0,
                     label='Choose the stage to analyze: ',
                     help='Choose the simulation stage to analyze'
                     )
      group.addParam('displayMdPymol', params.LabelParam,
                     label='Display trajectory with PyMol: ',
                     help='Display trajectory with Pymol. \n'
                          'Protein represented as NewCartoon and waters as sticks'
                     )
      group.addParam('displayMdVMD', params.LabelParam,
                     label='Display trajectory with VMD: ',
                     help='Display trajectory with VMD. \n'
                          'Protein represented as NewCartoon and waters as dots'
                     )

      group = form.addGroup('Gromacs analysis')
      group.addParam('displayAnalysis', params.EnumParam,
                     choices=self._analysis, default=0,
                     label='Choose the analysis to display: ',
                     help='Display the chosen analysis'
                     )
      line = group.addLine('Groups for analysis: ')
      line.addParam('chooseRef', params.EnumParam,
                     choices=self._ndxGroups, default=1,
                     label='Reference group: ',
                     help='Reference structure group to calculate the analysis against'
                     )
      line.addParam('chooseLsq', params.EnumParam,
                    choices=self._ndxGroups, default=1,
                    label='Least squares group: ', condition='displayAnalysis in [0, 5]',
                    help='Structure group to calculate the analysis in'
                    )
      group.addParam('aveRes', params.BooleanParam,
                    default=True,
                    label='Residue average: ', condition='displayAnalysis in [1]',
                    help='Display residue average'
                    )
      group.addParam('sasaOut', params.EnumParam,
                    choices=['Time', 'Residue', 'Atom'], default=0,
                    label='SASA average calculation over: ', condition='displayAnalysis in [3]',
                    help='Display sasa 1) group average over time 2) group average over residues '
                         '3) group average over atoms'
                    )
      line.addParam('chooseRef2', params.EnumParam,
                    choices=self._ndxGroups, default=0,
                    label='Reference group 2: ', condition='displayAnalysis in [4]',
                    help='Reference structure group 2 to calculate the analysis against'
                    )
      group.addParam('hbondOut', params.EnumParam,
                    choices=['Number', 'Autocorrelations', 'Distance', 'Angles'], default=0,
                    label='HBond option: ', condition='displayAnalysis in [4]',
                    help='Display Hbonds:'
                         '1) Number \n'
                         '2) Average over all autocorrelations of the existence functions (either 0 or 1) of all '
                         'hydrogen bonds. \n'
                         '3) Distances\n'
                         '4) Angles'
                    )
      group.addParam('clustMethod', params.EnumParam,
                     choices=['Single', 'Jarvis Patrick', 'Monte Carlo', 'Diagonalization', 'Gromos'],
                     default=4, label='Clustering method: ', condition='displayAnalysis in [5]',
                     help='Clustering method'
                     )
      group.addParam('clustCutoff', params.FloatParam,
                     default=0.1, label='Clustering cutoff: ', condition='displayAnalysis in [5]',
                     help='RMSD cut-off (nm) for two structures to be neighbor'
                     )

    def _getVisualizeDict(self):
      return {
        'displayPymol': self._showPymol,
        'displayMdPymol': self._showMdPymol,
        'displayMdVMD': self._showMdVMD,
        'displayAnalysis': self._showAnalysis,
      }

    def _showPymol(self, paramName=None):
      system = self.protocol.outputSystem

      return GromacsSystemViewer(project=self.getProject())._visualize(system)

    def _showMdPymol(self, paramName=None):
        stage = self.getEnumText('chooseStage')
        _, trjFile = self.getStageFiles(stage)
        system = self.protocol.outputSystem
        outPml = self.protocol._getExtraPath('pymolSimulation.pml')
        with open(outPml, 'w') as f:
          f.write(PML_MD_STR.format(os.path.abspath(system.getSystemFile()),
                                    os.path.abspath(trjFile)))

        return [PyMolView(os.path.abspath(outPml), cwd=self.protocol._getPath())]

    def _showMdVMD(self, paramName=None):
      stage = self.getEnumText('chooseStage')
      _, trjFile = self.getStageFiles(stage)
      system = self.protocol.outputSystem

      outTcl = self.protocol._getExtraPath('vmdSimulation.tcl')
      with open(outTcl, 'w') as f:
          f.write(TCL_MD_STR % (system.getSystemFile(), trjFile))
      args = '-e {}'.format(outTcl)

      return [VmdViewPopen(args)]

    def _showAnalysis(self, paramName=None):
        stage = self.getEnumText('chooseStage')
        if self.getEnumText('displayAnalysis') == 'RMSD':
            anFile = self.performRMSD(stage)
        elif self.getEnumText('displayAnalysis') == 'RMSF':
            anFile = self.performRMSF(stage)
        elif self.getEnumText('displayAnalysis') == 'Gyration':
            anFile = self.performGyration(stage)
        elif self.getEnumText('displayAnalysis') == 'SASA':
            anFile = self.performSASA(stage)
        elif self.getEnumText('displayAnalysis') == 'HBond':
            anFile = self.performHbond(stage)
        elif self.getEnumText('displayAnalysis') == 'Clustering':
            anFiles = self.performClustering(stage)

        if not self.getEnumText('displayAnalysis') in ['Clustering']:
            subprocess.Popen('xmgrace ' + anFile, shell=True)
        elif self.getEnumText('displayAnalysis') in ['Clustering']:
            print('Log file written in ', anFiles[1])
            modelsFiles = self.splitPDBModels(anFiles[0])
            modelSet = SetOfAtomStructs.create(self.getStageDir(stage))
            for mFile in modelsFiles:
                modelSet.append(AtomStruct(filename=mFile))

            chimViewer = ChimeraViewer(project=self.getProject(), protocol=self.protocol)
            return chimViewer._visualize(modelSet)


################################# ANALYSIS  #########################################

    def getCommandProgram(self):
      if self.getEnumText('displayAnalysis') in ['RMSD', 'Clustering']:
        return 'printf "{}\n{}\n" | '.format(self.getEnumText('chooseLsq'), self.getEnumText('chooseRef')) + program
      elif self.getEnumText('displayAnalysis') in ['RMSF', 'Gyration', 'SASA']:
        return 'printf "{}\n" | '.format(self.getEnumText('chooseRef')) + program
      elif self.getEnumText('displayAnalysis') in ['HBond']:
        return 'printf "{}\n{}\n" | '.format(self.getEnumText('chooseRef2'), self.getEnumText('chooseRef')) + program

    def performRMSD(self, stage):
        oFile = '%s_rmsd.xvg' % stage
        oDir = self.getStageDir(stage)
        oPath = os.path.join(oDir, oFile)
        if os.path.exists(oPath):
            os.remove(oPath)

        tprFile, trjFile = self.getStageFiles(stage)
        args = ' rms -s %s -f %s -o %s -tu ns' % (os.path.abspath(tprFile),
                                                  os.path.abspath(trjFile), oFile)
        subprocess.check_call(self.getCommandProgram() + args, shell=True, cwd=oDir)
        return oPath

    def performRMSF(self, stage):
        oFile = '%s_rmsf.xvg' % stage
        oDir = self.getStageDir(stage)
        oPath = os.path.join(oDir, oFile)
        if os.path.exists(oPath):
            os.remove(oPath)

        tprFile, trjFile = self.getStageFiles(stage)
        args = ' rmsf -s %s -f %s -o %s' % (os.path.abspath(tprFile), os.path.abspath(trjFile), oFile)
        if self.aveRes.get():
            args += ' -res'
        subprocess.check_call(self.getCommandProgram() + args, shell=True, cwd=oDir)
        return oPath

    def performGyration(self, stage):
        oFile = '%s_gyr.xvg' % stage
        oDir = self.getStageDir(stage)
        oPath = os.path.join(oDir, oFile)
        if os.path.exists(oPath):
            os.remove(oPath)

        tprFile, trjFile = self.getStageFiles(stage)
        args = ' gyrate -s %s -f %s -o %s' % (os.path.abspath(tprFile), os.path.abspath(trjFile), oFile)
        subprocess.check_call(self.getCommandProgram() + args, shell=True, cwd=oDir)
        return oPath

    def performSASA(self, stage):
        oFiles = ['%s_sasa.xvg' % stage, '%s_sasa_res.xvg' % stage, '%s_sasa_atom.xvg' % stage]
        outOptions = ['-o', '-or', '-oa']

        oFile = oFiles[self.sasaOut.get()]
        oDir = self.getStageDir(stage)
        oPath = os.path.join(oDir, oFile)
        if os.path.exists(oPath):
            os.remove(oPath)

        tprFile, trjFile = self.getStageFiles(stage)
        args = ' sasa -s %s -f %s %s %s -tu ns' % (os.path.abspath(tprFile), os.path.abspath(trjFile),
                                                   outOptions[self.sasaOut.get()], oFile)
        subprocess.check_call(self.getCommandProgram() + args, shell=True, cwd=oDir)
        return oPath

    def performHbond(self, stage):
        oFiles = ['%s_hbNum.xvg' % stage, '%s_hbAc.xvg' % stage, '%s_hbDist.xvg' % stage, '%s_hbAng.xvg' % stage]
        outOptions = ['-num', '-ac', '-dist', '-ang']

        oFile = oFiles[self.hbondOut.get()]
        oDir = self.getStageDir(stage)
        oPath = os.path.join(oDir, oFile)
        if os.path.exists(oPath):
            os.remove(oPath)

        tprFile, trjFile = self.getStageFiles(stage)
        args = ' hbond -s %s -f %s %s %s' % (os.path.abspath(tprFile), os.path.abspath(trjFile),
                                             outOptions[self.hbondOut.get()], oFile)
        subprocess.check_call(self.getCommandProgram() + args, shell=True, cwd=oDir)
        return oPath

    def performClustering(self, stage):
        oFiles = ['%s_cluster.pdb' % stage, '%s_cluster.log' % stage, 'rmsd-dist.xvg', 'rmsd-clust.xpm']
        oDir = self.getStageDir(stage)

        oPaths = []
        for oFile in oFiles:
            oPath = os.path.join(oDir, oFile)
            oPaths.append(oPath)
            if os.path.exists(oPath):
                os.remove(oPath)

        tprFile, trjFile = self.getStageFiles(stage)
        args = ' cluster -s %s -f %s -cl %s -g %s -method %s -cutoff %s' % \
               (os.path.abspath(tprFile), os.path.abspath(trjFile), oFiles[0], oFiles[1],
                self.getEnumText('clustMethod').lower(), self.clustCutoff.get())
        subprocess.check_call(self.getCommandProgram() + args, shell=True, cwd=oDir)
        return oPaths


################################# UTILS #########################

    def _getStagesWTrj(self):
        '''Return stages with a saved trajectory'''
        stages = ['All']
        for stDir in natural_sort(glob.glob(self.protocol._getExtraPath('stage_*'))):
            stage = os.path.basename(stDir)
            trjFile = '{}/{}.trr'.format(stDir, stage)
            if os.path.exists(trjFile):
                stages.append(stage)
        return stages

    def correctTrj(self, stage):
        comProg = 'printf "Protein\nSystem\n" | ' + program
        args = ' trjconv -s {}.tpr -f {}.trr -o {}_corrected.xtc -pbc mol -center'.format(*[stage]*3)
        subprocess.check_call(comProg + args, cwd=self.protocol._getExtraPath(stage), shell=True)
        return self.protocol._getExtraPath('{}/{}_corrected.xtc'.format(stage, stage))

    def getStageFiles(self, stage):
        if stage == 'All':
            _, _, tprFile = self.protocol.getPrevFinishedStageFiles(reverse=True)
            trjFile = self.protocol._getPath('outputTrajectory.xtc')
        else:
            _, _, tprFile = self.protocol.getPrevFinishedStageFiles(stage)
            trjFile = self.protocol._getExtraPath('{}/{}_corrected.xtc'.format(stage, stage))
            if not os.path.exists(trjFile):
              trjFile = self.correctTrj(stage)
        return tprFile, trjFile

    def splitPDBModels(self, combinedPDBFile):
        outFiles = []
        baseFile = combinedPDBFile.replace('.pdb', '_model%s.pdb')
        toWrite, i = '', 1
        with open(combinedPDBFile) as fIn:
            for line in fIn:
                if line.startswith('ENDMDL'):
                  outFiles.append(baseFile % i)
                  with open(outFiles[-1], 'w') as f:
                      f.write(toWrite)
                  toWrite = ''
                  i += 1

                else:
                    toWrite += line + '\n'
        return outFiles

    def getStageDir(self, stage):
      if stage == 'All':
        oDir = self.protocol._getExtraPath()
      else:
        oDir = self.protocol._getExtraPath(stage)
      return oDir





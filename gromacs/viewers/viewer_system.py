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

import os, glob
import pyworkflow.protocol.params as params

from pwem.objects import SetOfAtomStructs, AtomStruct
from pwem.viewers import ChimeraViewer, EmPlotter

from pwchem.viewers import VmdViewPopen, MDSystemViewer, MDSystemPViewer
from pwchem.utils import natural_sort
from pwchem.constants import TCL_MD_STR

from gromacs import Plugin as gromacsPlugin
from ..objects import GromacsSystem
from ..protocols import GromacsMDSimulation

program = gromacsPlugin.getGromacsBin()

class GromacsSystemPViewer(MDSystemPViewer):
    """ Visualize the output of Gromacs simulation """
    _label = 'Viewer Gromacs System'
    _targets = [GromacsSystem]

    def __init__(self, **args):
      super().__init__(**args)

    def getMDSystem(self, objType=GromacsSystem):
        if type(self.protocol) == objType:
            return self.protocol
        else:
            return self.protocol.outputSystem


class GromacsSimulationViewer(GromacsSystemPViewer):
    """ Visualize the output of Gromacs simulation """
    _label = 'Viewer Gromacs Simulation'
    _targets = [GromacsMDSimulation]
    _analysis = ['RMSD', 'RMSF', 'Gyration', 'SASA', 'HBond', 'Clustering']

    def __init__(self, **args):
      super().__init__(**args)

    def _defineParams(self, form):
        super()._defineParams(form)
        self._defineAnalysisParams(form)

    def _defineSimParams(self, form):
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
                            'Protein represented as NewCartoon and waters as dots')

    def _defineAnalysisParams(self, form):
      group = form.addGroup('Gromacs analysis')
      group.addParam('displayAnalysis', params.EnumParam,
                     choices=self._analysis, default=0,
                     label='Choose the analysis to display: ',
                     help='Display the chosen analysis'
                     )
      group.addParam('chain_name', params.EnumParam,
                     choices=self.getChainChoices(), default=0, condition='displayAnalysis in [1, 3]',
                     label='*Chain* to display analysis on: ',
                     help='Display the chosen analysis only in this chain'
                     )
      group.addParam('chooseRef', params.EnumParam,
                    choices=self.getIndexGroups(), default=1,
                    label='Reference group: ',
                    help='Reference structure group to calculate the analysis against'
                    )
      group.addParam('chooseLsq', params.EnumParam,
                    choices=self.getIndexGroups(), default=1,
                    label='Least squares group: ', condition='displayAnalysis in [0, 5]',
                    help='Structure group to calculate the analysis in'
                    )
      group.addParam('chooseRef2', params.EnumParam,
                     choices=self.getIndexGroups(), default=1,
                     label='Reference group 2: ', condition='displayAnalysis in [4]',
                     help='Reference structure group 2 to calculate the analysis against'
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
      return group

    def _getVisualizeDict(self):
        vDic = super()._getVisualizeDict()
        vDic.update({'displayAnalysis': self._showAnalysis})
        return vDic

    def _showMdPymol(self, paramName=None):
      stage = self.getEnumText('chooseStage')
      _, trjFile = self.getStageFiles(stage)
      system = self.getMDSystem()
      return MDSystemViewer(project=self.getProject())._visualize(system, trjFile=trjFile)

    def _showMdVMD(self, paramName=None):
      stage = self.getEnumText('chooseStage')
      _, trjFile = self.getStageFiles(stage)
      system = self.getMDSystem()

      outTcl = self.protocol._getExtraPath('vmdSimulation.tcl')
      sysExt = os.path.splitext(system.getSystemFile())[1][1:]
      trjExt = os.path.splitext(trjFile)[1][1:]
      self.writeTCL(outTcl, system.getSystemFile(), sysExt, trjFile, trjExt, system.getLigandTopologyFile())

      args = '-e {}'.format(outTcl)
      return [VmdViewPopen(args)]

    def _showAnalysis(self, paramName=None, saveFn=None):
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
        xs, ys, prevX = [[]], [[]], 0
        with open(anFile) as f:
          for line in f:
              if line.startswith('@'):
                  if line.split()[1] == 'title':
                      title = line.split('"')[-2]
                  elif line.split()[1] == 'xaxis':
                      xlabel = line.split('"')[-2]
                  elif line.split()[1] == 'yaxis':
                      ylabel = line.split('"')[-2]
              elif not line.startswith('#'):
                  xi = self.str2num(line.split()[0])
                  if prevX > xi:
                      xs.append([]), ys.append([])
                  xs[-1].append(xi)
                  ys[-1].append(float(line.split()[1]))
                  prevX = xi

        self.plotter = EmPlotter(x=1, y=1, windowTitle='Gromacs trajectory analysis')
        self.plotter.createSubPlot(title, xlabel, ylabel)
        if len(xs) > 1:
            system = self.getMDSystem()
            chainNames = system.getChainNames()
            for xsi, ysi, cn in zip(xs, ys, chainNames):
                if self.getEnumText('chain_name') in ['All', cn]:
                    self.plotter.plotData(xsi, ysi, '-', label=cn)
        else:
            self.plotter.plotData(xs[0], ys[0], '-')
        self.plotter.show()
        self.plotter.legend()
        if saveFn:
            self.plotter.savefig(saveFn)

      elif self.getEnumText('displayAnalysis') in ['Clustering']:
        print('Log file written in ', anFiles[1])
        modelsFiles = self.splitPDBModels(anFiles[0])
        modelSet = SetOfAtomStructs.create(self.getStageDir(stage))
        for mFile in modelsFiles:
          modelSet.append(AtomStruct(filename=mFile))

        chimViewer = ChimeraViewer(project=self.getProject(), protocol=self.protocol)
        return chimViewer._visualize(modelSet)

    ################################# ANALYSIS  #########################################

    def getIndexNDX(self, key):
      if key in ['RMSD', 'Clustering']:
        return [self.chooseLsq.get(), self.chooseRef.get()]
      elif key in ['RMSF', 'Gyration', 'SASA']:
        return [self.chooseRef.get()]
      elif key in ['HBond']:
        return [self.chooseRef2.get(), self.chooseRef.get()]

    def performRMSD(self, stage):
      oFile = '%s_rmsd.xvg' % stage
      oDir = self.getStageDir(stage)
      oPath = os.path.join(oDir, oFile)
      if os.path.exists(oPath):
        os.remove(oPath)

      groFile, trjFile = self.getStageFiles(stage)
      args = ' rms -s %s -f %s -o %s -tu ns' % (os.path.abspath(groFile), os.path.abspath(trjFile), oFile)
      if self.getIndexFile():
          args += f' -n {self.getIndexFile()}'
      gromacsPlugin.runGromacsPrintf(printfValues=self.getIndexNDX('RMSD'),
                                     args=args, cwd=oDir)
      return oPath

    def performRMSF(self, stage):
      oFile = '%s_rmsf.xvg' % stage
      oDir = self.getStageDir(stage)
      oPath = os.path.join(oDir, oFile)
      if os.path.exists(oPath):
        os.remove(oPath)

      groFile, trjFile = self.getStageFiles(stage)
      args = ' rmsf -s %s -f %s -o %s' % (os.path.abspath(groFile), os.path.abspath(trjFile), oFile)
      if self.getIndexFile():
          args += f' -n {self.getIndexFile()}'
      if self.aveRes.get():
        args += ' -res'

      gromacsPlugin.runGromacsPrintf(printfValues=self.getIndexNDX('RMSF'),
                                     args=args, cwd=oDir)
      return oPath

    def performGyration(self, stage):
      oFile = '%s_gyr.xvg' % stage
      oDir = self.getStageDir(stage)
      oPath = os.path.join(oDir, oFile)
      if os.path.exists(oPath):
        os.remove(oPath)

      groFile, trjFile = self.getStageFiles(stage)
      args = ' gyrate -s %s -f %s -o %s' % (os.path.abspath(groFile), os.path.abspath(trjFile), oFile)
      if self.getIndexFile():
          args += f' -n {self.getIndexFile()}'
      gromacsPlugin.runGromacsPrintf(printfValues=self.getIndexNDX('Gyration'),
                                     args=args, cwd=oDir)
      return oPath

    def performSASA(self, stage):
      oFiles = ['%s_sasa.xvg' % stage, '%s_sasa_res.xvg' % stage, '%s_sasa_atom.xvg' % stage]
      outOptions = ['-o', '-or', '-oa']

      oFile = oFiles[self.sasaOut.get()]
      oDir = self.getStageDir(stage)
      oPath = os.path.join(oDir, oFile)
      if os.path.exists(oPath):
        os.remove(oPath)

      groFile, trjFile = self.getStageFiles(stage)
      args = ' sasa -s %s -f %s %s %s -tu ns' % (os.path.abspath(groFile), os.path.abspath(trjFile),
                                                 outOptions[self.sasaOut.get()], oFile)
      if self.getIndexFile():
          args += f' -n {self.getIndexFile()}'
      gromacsPlugin.runGromacsPrintf(printfValues=self.getIndexNDX('SASA'),
                                     args=args, cwd=oDir)
      return oPath

    def performHbond(self, stage):
      oFiles = ['%s_hbNum.xvg' % stage, '%s_hbAc.xvg' % stage, '%s_hbDist.xvg' % stage, '%s_hbAng.xvg' % stage]
      outOptions = ['-num', '-ac', '-dist', '-ang']

      oFile = oFiles[self.hbondOut.get()]
      oDir = self.getStageDir(stage)
      oPath = os.path.join(oDir, oFile)
      if os.path.exists(oPath):
        os.remove(oPath)

      tprFile, trjFile = self.getStageFiles(stage, tpr=True)
      args = ' hbond -s %s -f %s %s %s' % (os.path.abspath(tprFile), os.path.abspath(trjFile),
                                           outOptions[self.hbondOut.get()], oFile)
      if self.getIndexFile():
          args += f' -n {self.getIndexFile()}'
      gromacsPlugin.runGromacsPrintf(printfValues=self.getIndexNDX('HBond'),
                                     args=args, cwd=oDir)
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

      groFile, trjFile = self.getStageFiles(stage)
      args = ' cluster -s %s -f %s -cl %s -g %s -method %s -cutoff %s' % \
             (os.path.abspath(groFile), os.path.abspath(trjFile), oFiles[0], oFiles[1],
              self.getEnumText('clustMethod').lower(), self.clustCutoff.get())
      if self.getIndexFile():
          args += f' -n {self.getIndexFile()}'
      gromacsPlugin.runGromacsPrintf(printfValues=self.getIndexNDX('Clustering'),
                                     args=args, cwd=oDir)
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
      args = ' trjconv -s {}.tpr -f {}.trr -o {}_corrected.xtc -pbc mol -center'.format(*[stage] * 3)
      gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'System'],
                                     args=args, cwd=self.protocol._getExtraPath(stage))
      return self.protocol._getExtraPath('{}/{}_corrected.xtc'.format(stage, stage))

    def getStageFiles(self, stage, tpr=False):
      if stage == 'All':
        system = self.getMDSystem()
        groFile, trjFile, tprFile = system.getOriStructFile(), system.getTrajectoryFile(), system.getTprFile()
      else:
        groFile, _, tprFile = self.protocol.getPrevFinishedStageFiles(stage)
        trjFile = self.protocol._getExtraPath('{}/{}_corrected.xtc'.format(stage, stage))
        if not os.path.exists(trjFile):
          trjFile = self.correctTrj(stage)
      if not tpr:
          return groFile, trjFile
      else:
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

    def str2num(self, stri):
        x = float(stri)
        try:
            x2 = int(x)
            if x2 == x:
                return x2
            else:
                return x
        except:
            pass

    def getChainChoices(self):
        system = self.getMDSystem()
        return ['All'] + system.getChainNames()

    def getIndexFile(self):
        indexFile = self.getMDSystem().getIndexFile()
        if indexFile and os.path.exists(indexFile):
            return os.path.abspath(indexFile)

    def getIndexGroupsDic(self):
        groups = self.protocol.parseIndexFile(self.protocol.getCustomIndexFile())
        return groups

    def getIndexGroups(self):
        groups = self.getIndexGroupsDic()
        return list(groups.values())


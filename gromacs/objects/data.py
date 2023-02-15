# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Pedro Febrer Martinez (pedrofebrer98@gmail.com)
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

import os, shutil
from subprocess import check_call
import pyworkflow.object as pwobj

from pwchem.objects import MDSystem
from gromacs.constants import *

class GromacsSystem(MDSystem):
    """A system atom structure (prepared for MD) in the file format of GROMACS
    _topoFile: topology file .top
    _restrFile: default position restraints file .itp
    _trjFile: trajectory file (xtc)
    _ff: main force field
    _wff: water force field model"""

    def __init__(self, filename=None, **kwargs):
        super().__init__(filename=filename, **kwargs)
        self._oriStructFile = pwobj.String(kwargs.get('oriStructFile', None))
        self._restrFile = pwobj.String(kwargs.get('restrFile', None))
        self._tprFile = pwobj.String(kwargs.get('tprFile', None))

        self._chainNames = pwobj.String(kwargs.get('chainNames', None))

        self._firstFrame = pwobj.Integer(kwargs.get('firstFrame', None))
        self._lastFrame = pwobj.Integer(kwargs.get('lastFrame', None))

        self._firstTime = pwobj.Float(kwargs.get('firstTime', None))
        self._lastTime = pwobj.Float(kwargs.get('lastTime', None))

    def __str__(self):
        strStr = '{} ({}'.format(self.getClassName(), os.path.basename(self.getSystemFile()))
        if self.hasTrajectory():
            strStr += ', frames: {} - {}, time(ps): {} - {}'.format(*self.getFrameIdxs(), *self.getTimes())
        strStr += ')'
        return strStr

    def getChainNames(self):
        return self._chainNames.get().split(',')
    def setChainNames(self, values):
        if type(values) == str:
            self._chainNames.set(values)
        elif type(values) in [list, tuple]:
            self._chainNames.set(','.join(values))

    def getFrameIdxs(self):
        return self._firstFrame, self._lastFrame
    def setFrameIdxs(self, values):
        self._firstFrame.set(values[0]), self._lastFrame.set(values[1])

    def getTimes(self):
        return self._firstTime, self._lastTime
    def setTimes(self, values):
        self._firstTime.set(values[0]), self._lastTime.set(values[1])

    def readTrjInfo(self, protocol, outDir=None):
        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir
        infoFile = os.path.abspath(protocol._getPath('logs/run.stderr'))

        command = 'check -f {}'.format(self.getTrajectoryFile())
        gromacsPlugin.runGromacs(protocol, 'gmx', command, cwd=outDir)

        isCheck, isFirst = True, True
        with open(infoFile) as f:
          for line in f:
              if isCheck:
                  if 'gmx check -f' in line:
                      isCheck = False
              else:
                  if isFirst and line.startswith('Reading frame'):
                      isFirst = False
                      firstFrame, firstTime = int(line.split()[2]), float(line.split()[4])
                  elif not isFirst and line.startswith('Last frame'):
                      lastFrame, lastTime = int(line.split()[2]), float(line.split()[4])

        self.setTimes([firstTime, lastTime])
        self.setFrameIdxs([firstFrame, lastFrame])

    def getOriStructFile(self):
        return self._oriStructFile.get()

    def setOriStructFile(self, value):
        self._oriStructFile.set(value)

    def getRestraintsFile(self):
        return self._restrFile.get()

    def setRestraintsFile(self, value):
        self._restrFile.set(value)

    def defineNewRestrictionWrapper(self, energy, restrainSuffix='low', outDir=None, group="protein-h"):
        '''Call the right function to define a new position restriction and 
        store it in the right topology files, checking whether there is one chain
        or more than one.'''

        # check if there are topol itp files
        lastTopDir = os.path.dirname(self.getTopologyFile())
        lastTopItpFiles = [os.path.join(lastTopDir, f)
                            for f in os.listdir(lastTopDir)
                            if f.startswith('top') and f.endswith('.itp')]

        if len(lastTopItpFiles) == 0:
            return [self.defineNewRestriction(energy, restrainSuffix, outDir, group)]
        else:
            return self.defineNewRestrictionMulti(energy, restrainSuffix, outDir, 
                                                  group, lastTopItpFiles)

    def defineNewRestrictionMulti(self, energy, restrainSuffix='low', outDir=None, 
                                  group="protein-h", lastTopItpFiles=[]):
        '''Define a new position restriction and store it in the right topology files
        for each of the multiple chains.'''
        lastTopItpFiles.sort(key=os.path.getmtime)

        topFile = self.getTopologyFile()
        if outDir:
            topFile = os.path.join(outDir, os.path.basename(self.getTopologyFile()))
            shutil.copy(self.getTopologyFile(), topFile)
            self.setTopologyFile(topFile)
        
        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir

        if group == "protein-h":
            groupNr = 2
        else:
            print("Using protein-h as other groups not supported")
            groupNr = 2

        # make standard ndx file with simple quitting
        paramsFileName = os.path.abspath(os.path.join(outDir, "make_ndx_inputs.txt"))
        fo = open(paramsFileName, "w")
        fo.write("q\n")
        fo.close()

        program = os.path.join("", '{} '.format(gromacsPlugin.getGromacsBin()))
        ndx_file = os.path.abspath(os.path.join(outDir, 'index.ndx'))
        params_make_ndx = 'make_ndx -f %s -o %s < %s' % \
                        (os.path.abspath(self.getSystemFile()), 
                        ndx_file, paramsFileName)
        check_call(program + params_make_ndx, cwd=outDir, shell=True)

        # count ndx entries
        f = open(ndx_file, 'r')
        numGroups = 0
        for line in f.readlines():
            if line.find('[') != -1:
                numGroups += 1
        f.close()

        # make ndx file with split chains from group
        paramsFileName = os.path.abspath(os.path.join(outDir, "make_ndx_inputs.txt"))
        fo = open(paramsFileName, "w")
        fo.write("splitch {}\nq\n".format(groupNr))
        fo.close()

        program = os.path.join("", '{} '.format(gromacsPlugin.getGromacsBin()))
        ndx_file = os.path.abspath(os.path.join(outDir, 'index.ndx'))
        params_make_ndx = 'make_ndx -f %s -n %s -o %s < %s' % \
                        (os.path.abspath(self.getSystemFile()), 
                         ndx_file, ndx_file, paramsFileName)
        check_call(program + params_make_ndx, cwd=outDir, shell=True)

        suffixes = []
        for n, TopItpFile in enumerate(lastTopItpFiles):
            # copy the file
            itpFile = os.path.join(outDir, os.path.basename(TopItpFile))
            shutil.copy(TopItpFile, itpFile)

            # generate restraint for the new group for this chain
            newRestrainSuffix = '%s_%d' % (restrainSuffix.lower(), numGroups+n)
            program = os.path.join("", 'printf "{}" | {} '.format(numGroups+n, gromacsPlugin.getGromacsBin()))
            params_genrestr = 'genrestr -f %s -o posre_%s.itp -fc %d %d %d -n %s' % \
                            (os.path.abspath(self.getSystemFile()),
                            newRestrainSuffix, 
                            energy, energy, energy, ndx_file)
            check_call(program + params_genrestr, cwd=outDir, shell=True)  

            program = "sed "
            inStr = '#ifdef POSRES_{}\\n#include "posre_{}.itp"\\n#endif'.\
                format(newRestrainSuffix.upper(), newRestrainSuffix.lower())
            sed_params = """-i '/; Include Position restraint file/a {}' {}""".\
                format(inStr, os.path.abspath(itpFile))
            check_call(program + sed_params, cwd=outDir, shell=True)   

            suffixes.append(newRestrainSuffix)

        return suffixes

    def defineNewRestriction(self, energy, restrainSuffix='low', outDir=None, group="protein-h"):
        '''Define a new position restriction and store it in the topology file.
        This only works for simple systems with one chain'''

        if group == "protein-h":
            groupNr = 2
        else:
            print("Using protein-h as other groups not supported")
            groupNr = 2

        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir
        program = os.path.join("", 'printf "{}" | {} '.format(groupNr, gromacsPlugin.getGromacsBin()))
        params_genrestr = 'genrestr -f %s -o %s.itp -fc %d %d %d' % \
                          (os.path.abspath(self.getSystemFile()),
                           'posre_' + restraintSuffix.lower(), energy, energy, energy)
        check_call(program + params_genrestr, cwd=outDir, shell=True)

        topFile = self.getTopologyFile()
        if outDir:
            topFile = os.path.join(outDir, os.path.basename(self.getTopologyFile()))
            shutil.copy(self.getTopologyFile(), topFile)
            self.setTopologyFile(topFile)

        program = "sed "
        inStr = '#ifdef POSRES_{}\\n#include "posre_{}.itp"\\n#endif'.\
            format(restraintSuffix.upper(), restraintSuffix.lower())
        sed_params = """-i '/; Include Position restraint file/a {}' {}""".\
            format(inStr, os.path.abspath(topFile))
        check_call(program + sed_params, cwd=outDir, shell=True)

        return restrainSuffix

    def getIons(self):
        ionsDic, mols = {}, False
        with open(self.getTopologyFile()) as f:
            for line in f:
                if mols:
                    sline = line.split()
                    if sline[0] in ION_NAMES:
                        ionsDic[sline[0]] = int(sline[1])
                if line.startswith('[ molecules ]'):
                    mols = True
        return ionsDic



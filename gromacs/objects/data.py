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
        self._restrFile = pwobj.String(kwargs.get('restrFile', None))
        self._tprFile = pwobj.String(kwargs.get('tprFile', None))
        self._indexFile = pwobj.String(kwargs.get('indexFile', None))

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
        if isinstance(values, str):
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

        command = 'check -f {}'.format(os.path.abspath(self.getTrajectoryFile()))
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
                  elif not isFirst and line.startswith('Reading frame'):
                      lastFrame, lastTime = int(line.split()[2]), float(line.split()[4])

        self.setTimes([firstTime, lastTime])
        self.setFrameIdxs([firstFrame, lastFrame])

    def getRestraintsFile(self):
        return self._restrFile.get()

    def setRestraintsFile(self, value):
        self._restrFile.set(value)

    def getTprFile(self):
        return self._tprFile.get()

    def setTprFile(self, value):
        self._tprFile.set(value)

    def getIndexFile(self):
        return self._indexFile.get()

    def setIndexFile(self, value):
        self._indexFile.set(value)

    def defineNewRestriction(self, index, energy, restraintSuffix='low', outDir=None, indexFile=None):
        '''Define a new position restriction and stores it in the topology file'''
        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir

        nArg = ' -n {}'.format(indexFile) if indexFile else ''
        paramsGenrestr = 'genrestr -f %s%s -o %s.itp -fc %d %d %d' % \
                          (os.path.abspath(self.getSystemFile()), nArg,
                           'posre_' + restraintSuffix.lower(), energy, energy, energy)
        gromacsPlugin.runGromacsPrintf(printfValues=index, args=paramsGenrestr, cwd=outDir)

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



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
import pwem.objects.data as data

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
        self._oriStructFile = pwobj.String(kwargs.get('oriStructFile', None))

        self._freeEnergy = pwobj.Float(kwargs.get('MMPGSA', None))
        self._freeEnergyFile = pwobj.String(kwargs.get('MMPGSAFile', None))
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
        value = os.path.relpath(value)
        self._indexFile.set(value)

    def getFreeEnergy(self):
        return self._freeEnergy.get()

    def setFreeEnergy(self, value):
        value = float(value)
        self._freeEnergy.set(value)

    def getFreeEnergyFile(self):
        return self._freeEnergyFile.get()

    def setFreeEnergyFile(self, value):
        value = os.path.relpath(value)
        self._freeEnergyFile.set(value)

    def getOriStructFile(self):
        return self._oriStructFile.get()

    def setOriStructFile(self, value):
        self._oriStructFile.set(value)

    def defineNewRestriction(self, protocol, index, energy, restraintSuffix='low', outDir=None, indexFile=None):
        '''Define a new position restriction and stores it in the topology file'''
        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir

        nArg = ' -n {}'.format(indexFile) if indexFile else ''
        paramsGenrestr = 'genrestr -f %s%s -o %s.itp -fc %d %d %d' % \
                          (os.path.abspath(self.getSystemFile()), nArg,
                           'posre_' + restraintSuffix.lower(), energy, energy, energy)
        gromacsPlugin.runGromacsPrintf(protocol, printfValues=index, args=paramsGenrestr, cwd=outDir)

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


class FreeEnergyCalculation(data.EMFile):
    """
    Object representing the results of a Binding Free Energy calculation
    (such as MM/PBSA or MM/GBSA) for a protein-ligand complex.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        # Core values
        self._deltaG = pwobj.Float(kwargs.get('deltaG', None))
        self._deviation = pwobj.Float(kwargs.get('deviation', None))
        self._calcType = pwobj.String(kwargs.get('calcType', None))

        # Output files
        self._resultFile = pwobj.String(kwargs.get('resultFile', None))
        self._csvFile = pwobj.String(kwargs.get('csvFile', None))

        # Metadata
        self._entropyType = pwobj.String(kwargs.get('entropyType', 'None'))

    def __str__(self):
        calcType = self.getCalcType()
        dg = self.getDeltaG()
        dgStr = f"{dg:.2f} kcal/mol" if dg is not None else "N/A"
        return f"FreeEnergyCalculation ({calcType}  (ΔG = {dgStr})"

    def getDeltaG(self):
        return self._deltaG.get()

    def setDeltaG(self, value):
        self._deltaG.set(float(value))

    def getCalcType(self):
        return self._calcType.get()

    def setCalcType(self, value):
        self._calcType.set(value)

    def getResultFile(self):
        return self._resultFile.get()

    def setResultFile(self, value):
        if value is not None:
            value = os.path.relpath(value)
        self._resultFile.set(value)

    def getCsvFile(self):
        return self._csvFile.get()

    def setCsvFile(self, value):
        if value is not None:
            value = os.path.relpath(value)
        self._csvFile.set(value)

    def getEntropyType(self):
        return self._entropyType.get()

    def setEntropyType(self, value):
        self._entropyType.set(value)

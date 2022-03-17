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
import pwem.objects.data as data
import pyworkflow.object as pwobj

class GromacsSystem(data.EMFile):
    """A system atom structure (prepared for MD) in the file format of GROMACS
    _topoFile: topology file .top
    _restrFile: default position restrains file .itp
    _trjFile: trajectory file (xtc)
    _ff: main force field
    _wff: water force field model"""

    def __init__(self, filename=None, **kwargs):
        super().__init__(filename=filename, **kwargs)
        self._topoFile = pwobj.String(kwargs.get('topoFile', None))
        self._restrFile = pwobj.String(kwargs.get('restrFile', None))
        self._trjFile = pwobj.String(kwargs.get('trjFile', None))
        self._ff = pwobj.String(kwargs.get('ff', None))
        self._wff = pwobj.String(kwargs.get('wff', None))

    def __str__(self):
        return '{} ({}, hasTrj={})'.format(self.getClassName(), os.path.basename(self.getSystemFile()),
                                           self.hasTrajectory())

    def getSystemFile(self):
        return self.getFileName()

    def setSystemFile(self, value):
        self.setFileName(value)

    def getTopologyFile(self):
        return self._topoFile.get()

    def setTopologyFile(self, value):
        self._topoFile.set(value)

    def getRestrainsFile(self):
        return self._restrFile.get()

    def setRestrainsFile(self, value):
        self._restrFile.set(value)

    def hasTrajectory(self):
        if self.getTrajectoryFile():
            return True
        else:
            return False

    def getTrajectoryFile(self):
        return self._trjFile.get()

    def setTrajectoryFile(self, value):
        self._trjFile.set(value)

    def getForceField(self):
        return self._ff

    def setForceField(self, value):
        self._ff.set(value)

    def getWaterForceField(self):
        return self._wff

    def setWaterForceField(self, value):
        self._wff.set(value)

    def defineNewRestriction(self, energy, restrainSuffix='low', outDir=None):
        '''Define a new position restriction and stores it in the topology file'''
        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir
        program = os.path.join("", 'printf "2" | {} '.format(gromacsPlugin.getGromacsBin()))
        params_genrestr = 'genrestr -f %s -o %s.itp -fc %d %d %d' % \
                          (os.path.abspath(self.getSystemFile()),
                           'posre_' + restrainSuffix.lower(), energy, energy, energy)
        check_call(program + params_genrestr, cwd=outDir, shell=True)

        topFile = self.getTopologyFile()
        if outDir:
            topFile = os.path.join(outDir, os.path.basename(self.getTopologyFile()))
            shutil.copy(self.getTopologyFile(), topFile)
            self.setTopologyFile(topFile)

        program = "sed "
        inStr = '#ifdef POSRES_{}\\n#include "posre_{}.itp"\\n#endif'.\
            format(restrainSuffix.upper(), restrainSuffix.lower())
        sed_params = """-i '/; Include Position restraint file/a {}' {}""".\
            format(inStr, os.path.abspath(topFile))
        check_call(program + sed_params, cwd=outDir, shell=True)

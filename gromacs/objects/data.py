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

    def getRestraintsFile(self):
        return self._restrFile.get()

    def setRestraintsFile(self, value):
        self._restrFile.set(value)

    def defineNewRestriction(self, index, energy, restraintSuffix='low', outDir=None):
        '''Define a new position restriction and stores it in the topology file'''
        from gromacs import Plugin as gromacsPlugin
        outDir = os.path.dirname(self.getSystemFile()) if not outDir else outDir
        program = os.path.join("", 'printf "{}" | {} '.format(index, gromacsPlugin.getGromacsBin()))
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



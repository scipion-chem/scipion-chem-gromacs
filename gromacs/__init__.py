# **************************************************************************
# *
# * Authors:      Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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
import pwem
from .constants import *
from os.path import join
import subprocess

_logo = "icon.png"
_references = ['Abraham2015']

from .objects import *
class Plugin(pwem.Plugin):
    _homeVar = GROMACS_HOME
    _pathVars = [GROMACS_HOME]
    _supportedVersions = [V2020]
    _gromacsName = GROMACS + '-' + GROMACS_DEFAULT_VERSION
    _pluginHome = join(pwem.Config.EM_ROOT, _gromacsName)

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(GROMACS_HOME, cls._gromacsName)

    @classmethod
    def defineBinaries(cls, env):
        addC36cmd = 'cd share/top && wget -O charmm36-feb2021.ff.tgz http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-feb2021.ff.tgz && '
        addC36cmd += 'tar -xf charmm36-feb2021.ff.tgz'
        cMakeCmd = 'mkdir build && cd build && '
        cMakeCmd += 'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA ' \
                    '-DCMAKE_INSTALL_PREFIX={}  -DGMX_FFT_LIBRARY=fftw3 > cMake.log'.format(cls._pluginHome)
        makeCmd = 'cd build && make -j {} > make.log && make check > check.log'.format(env.getProcessors())
        makeInstallCmd = 'cd build && make install > install.log'

        # Creating validation file
        GROMACS_INSTALLED = '%s_installed' % GROMACS
        installationCmd = 'touch %s' % GROMACS_INSTALLED  # Flag installation finished

        env.addPackage(GROMACS,
                       version=GROMACS_DEFAULT_VERSION,
                       url=cls._getGromacsDownloadUrl(),
                       commands=[(addC36cmd, 'share/top/charmm36-feb2021.ff'), 
                                 (cMakeCmd, 'build/cMake.log'),
                                 (makeCmd, 'build/check.log'),
                                 (makeInstallCmd, 'build/install.log'),
                                 (installationCmd, GROMACS_INSTALLED)],
                       default=True)

    @classmethod
    def runGromacs(cls, protocol, program, args, cwd=None):
        """ Run Gromacs command from a given protocol. """
        protocol.runJob(cls.getGromacsBin(program), args, cwd=cwd)

    @classmethod
    def runGromacsPrintf(cls, printfValues, args, cwd):
      """ Run Gromacs command from a given protocol. """
      program = 'printf "{}\n" | {} '.format('\n'.join(printfValues), cls.getGromacsBin())
      subprocess.check_call(program + args, cwd=cwd, shell=True)

    @classmethod
    def getGromacsBin(cls, program='gmx'):
        return join(cls._pluginHome, 'bin/{}'.format(program))

    @classmethod  # Test that
    def getEnviron(cls):
        pass

    @classmethod
    def _getGromacsDownloadUrl(cls):
        return 'https://ftp.gromacs.org/gromacs/gromacs-{}.tar.gz'.format(GROMACS_DEFAULT_VERSION)

    @classmethod
    def _getGromacsTar(cls):
        return cls._pluginHome + '/' + cls._gromacsName + '.tar.gz'

    # ---------------------------------- Utils functions  -----------------------


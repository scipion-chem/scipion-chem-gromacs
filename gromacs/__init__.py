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

_logo = "gromacs_logo.png"
_references = ['Abraham2015']
V2020 = '2020.6'
V2021 = '2021.5'

GROMACS_DIC = {'name': 'gromacs', 'version': '2021.5', 'home': 'GROMACS_HOME'}

from .objects import *
class Plugin(pwem.Plugin):
    _homeVar = GROMACS_DIC['home']
    _pathVars = [GROMACS_DIC['home']]
    _supportedVersions = [V2020, V2021]
    _gromacsName = GROMACS_DIC['name'] + '-' + GROMACS_DIC['version']

    @classmethod
    def _defineVariables(cls):
        """ Return and write a variable in the config file.
        """
        cls._defineEmVar(GROMACS_DIC['home'], cls._gromacsName)

    @classmethod
    def defineBinaries(cls, env):
        installationCmd = 'wget -O gromacs-{}.tar.gz --no-check-certificate {} && '. \
          format(GROMACS_DIC['version'], cls._getGromacsDownloadUrl())
        installationCmd += 'tar -xf gromacs-{}.tar.gz --strip-components 1 && '.\
          format(GROMACS_DIC['version'], GROMACS_DIC['version'])
        installationCmd += 'cd share/top && wget -O charmm36-feb2021.ff.tgz http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-feb2021.ff.tgz && '
        installationCmd += 'tar -xf charmm36-feb2021.ff.tgz && cd ../.. && '
        installationCmd += 'mkdir build && cd build && '

        cmakeCmd = 'cmake'
        import subprocess
        cmake_version = subprocess.check_output(['cmake', '--version']).decode().split('\n')[0].split()[-1]
        if cmake_version.startswith('3'):
            cmakeCmd = 'cmake'
        else:
            try:
                cmake_version = subprocess.check_output(['cmake3', '--version']).decode().split('\n')[0].split()[-1]
                if cmake_version.startswith('3') and int(cmake_version.split('.')[1]) > 13:
                    cmakeCmd = 'cmake3'
                else:
                    print('cmake3 is not higher than 3.13')
            except:
                print('cmake is not higher than 3.13')

        installationCmd += cmakeCmd + ' .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA ' \
                    '-DCMAKE_INSTALL_PREFIX={}  -DGMX_FFT_LIBRARY=fftw3 > cMake.log && '.\
          format(cls.getVar(GROMACS_DIC['home']))
        installationCmd += 'make -j {} > make.log && make check > check.log && '.format(env.getProcessors())
        installationCmd += 'make install > install.log && '

        # Creating validation file
        GROMACS_INSTALLED = '%s_installed' % GROMACS_DIC['name']
        installationCmd += 'cd .. && touch %s' % GROMACS_INSTALLED  # Flag installation finished

        env.addPackage(GROMACS_DIC['name'],
                       version=GROMACS_DIC['version'],
                       tar='void.tgz',
                       commands=[(installationCmd, GROMACS_INSTALLED)],
                       default=True)

    @classmethod
    def runGromacs(cls, protocol, program='gmx', args='', cwd=None):
        """ Run Gromacs command from a given protocol. """
        protocol.runJob(cls.getGromacsBin(program), args, cwd=cwd)

    @classmethod
    def runGromacsPrintf(cls, printfValues, args, cwd):
      """ Run Gromacs command from a given protocol. """
      program = 'printf "{}\n" | {} '.format('\n'.join(printfValues), cls.getGromacsBin())
      print('Running: ', program, args)
      subprocess.check_call(program + args, cwd=cwd, shell=True)

    @classmethod
    def getGromacsBin(cls, program='gmx'):
        return join(cls.getVar(GROMACS_DIC['home']), 'bin/{}'.format(program))

    @classmethod  # Test that
    def getEnviron(cls):
        pass

    @classmethod
    def _getGromacsDownloadUrl(cls):
        return 'https://ftp.gromacs.org/gromacs/gromacs-{}.tar.gz'.format(GROMACS_DIC['version'])

    # ---------------------------------- Utils functions  -----------------------


# **************************************************************************
# *
# * Authors:	  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

# General imports
import os
import subprocess, multiprocessing
import os
from os.path import join

# Scipion em imports
import pwem
from scipion.install.funcs import InstallHelper
from pyworkflow.utils import redStr, yellowStr, Environ

# Plugin imports
from pwchem import Plugin as pwchemPlugin

from .objects import *
from .constants import *

_logo = "gromacs_logo.png"
__version__ = '1.0.0'
_references = ['Abraham2015']

class Plugin(pwchemPlugin):
	_homeVar = GROMACS_DIC['home']
	_pathVars = [GROMACS_DIC['home']]
	_supportedVersions = [V2020, V2021, V2026]
	_gromacsName = GROMACS_DIC['name'] + '-' + GROMACS_DIC['version']

	@classmethod
	def _defineVariables(cls):
		""" Return and write a variable in the config file. """
		cls._defineEmVar(GROMACS_DIC['home'], cls._gromacsName)

	@classmethod
	def defineBinaries(cls, env):
		""" This function defines all the packages that will be installed. """
		# Checking requirements
		modifiedProcs = cls.checkRequirements(env)

		# Installing packages
		cls.addGromacs(env, modifiedProcs)

	@classmethod
	def addGromacs(cls, env, modifiedProcs, default=True):
		""" This function installs Gromacs's package. """
		# Target suffix for MPI installation
		mpiExt = '_MPI'

		# Instantiating install helper
		installer = InstallHelper(GROMACS_DIC['name'], cls.getVar(GROMACS_DIC['home']), GROMACS_DIC['version'])

		# Defining some variables
		gromacsFileName = f"gromacs-{GROMACS_DIC['version']}.tar.gz"
		
		charmInnerLocation = join('share', 'top')
		charmFileName = 'charmm36-feb2021.ff.tgz'

		normalInnerLocation = 'build'
		mpiInnerLocation = 'build_mpi'

		# If number of processors has been modified, show message
		if modifiedProcs:
			message = "\\nWARNING: Only 1 process has been defined to install Gromacs.\\n"
			message += f"This will take a very long time. Instead, the number of parallel processes has been changed to the maximum avaliable in your system: {env.getProcessors()}."
			installer.addCommand(f'echo -e "{yellowStr(message)}"', 'WARNING_SHOWN')

		# Installing package
		installer.getExtraFile(cls._getGromacsDownloadUrl(), 'GROMACS_DOWNLOADED', fileName=gromacsFileName)\
			.addCommand(f'tar -xf {gromacsFileName} --strip-components 1', 'GROMACS_EXTRACTED')\
			.getExtraFile('http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-feb2021.ff.tgz', 'CHARM_DOWNLOADED', location=charmInnerLocation, fileName=charmFileName)\
			.addCommand(f'tar -xf {charmFileName}', 'CHARM_EXTRACTED', workDir=charmInnerLocation)\
			.addCommand(f'mkdir {normalInnerLocation} {mpiInnerLocation}', 'BUILD_DIRS_MADE')\
			.addCommand(f'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_CUDA_ARCHITECTURES=native -DCMAKE_INSTALL_PREFIX={cls.getVar(GROMACS_DIC["home"])}/install -DGMX_FFT_LIBRARY=fftw3',
	       				'GROMACS_BUILT', workDir=normalInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()}', 'GROMACS_COMPILED', workDir=normalInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()} install', 'GROMACS_INSTALLED', workDir=normalInnerLocation)\
			.addCommand(f'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX={cls.getVar(GROMACS_DIC["home"])}/install{mpiExt.lower()} -DGMX_FFT_LIBRARY=fftw3 -DGMX_MPI=on', 
	       				'GROMACS_BUILT' + mpiExt, workDir=mpiInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()}', 'GROMACS_COMPILED' + mpiExt, workDir=mpiInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()} install', 'GROMACS_INSTALLED' + mpiExt, workDir=mpiInnerLocation)\
			.addPackage(env, dependencies=['wget', 'tar', 'cmake', 'make'], default=default)

	@classmethod
	def runGromacs(cls, protocol, program='gmx', args='', cwd=None, mpi=False, **kwargs):
		""" Run Gromacs command from a given protocol. """
		protocol.runJob(cls.getGromacsBin(program, mpi=mpi), args, cwd=cwd, **kwargs)

	@classmethod
	def runGromacsPrintf(cls, protocol, printfValues, args, cwd, mpi=False):
		""" Run Gromacs command with interactive printf input via Scipion's runJob. """
		printfValues = list(map(str, printfValues))
		gmxBin = cls.getGromacsBin(mpi=mpi)
		fullProgram = 'printf "{}\\n" | {}'.format('\\n'.join(printfValues), gmxBin)

		protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd,
		                numberOfMpi=1, numberOfThreads=1)

	@classmethod
	def runGromacsPrintfViewer(cls, printfValues, args, cwd, mpi=False):
		""" Run Gromacs command from a given protocol. """
		printfValues = list(map(str, printfValues))
		program = 'printf "{}\n" | {} '.format('\n'.join(printfValues), cls.getGromacsBin(mpi=mpi))
		print('Running: ', program, args)
		subprocess.check_call(program + args, cwd=cwd, shell=True)

	@classmethod
	def getGromacsBin(cls, program='gmx', mpi=False):
		mpiExt = '_mpi' if mpi else ''
		return join(cls.getVar(GROMACS_DIC['home']), f'install{mpiExt}/bin/{program}{mpiExt}')

	@classmethod
	def getEnviron(cls):
		return Environ(os.environ)

	@classmethod
	def _getGromacsDownloadUrl(cls):
		return 'https://ftp.gromacs.org/gromacs/gromacs-{}.tar.gz'.format(GROMACS_DIC['version'])

	# ---------------------------------- Utils functions  -----------------------
	@classmethod
	def checkRequirements(cls, env):
		""" This function checks if the software requirements are being met. """
		cls.checkCMakeVersion()
		cls.checkCudaVersion()
		return cls.defineProcessors(env)

	@classmethod
	def defineProcessors(cls, env):
		""" This function defines the number of processors that will be used during installation. """
		# Check if user defined number of processes is 1
		# If so, set that number to the number of processes available
		if env.getProcessors() == 1:
			env._processors = multiprocessing.cpu_count()
			return True

	@classmethod
	def versionTuple(cls, versionStr):
		"""
		This function returns the given version sting ('1.0.7' for example) into a tuple, so it can be compared.
		It also accepts other version schemes, like 1.0.9-rc, but only the numeric part is taken into account.
		"""
		# Split the version string by dots
		versionParts = versionStr.split('.')
		# Initialize an empty list to store the numerical parts of the version string
		numericalParts = []
		# Iterate over each part of the version string
		for part in versionParts:
				# Split the part by hyphens
				subparts = part.split('-')
				# The first subpart is always numerical, so we append it to our list
				numericalParts.append(int(subparts[0]))
		# Convert the list of numerical parts to a tuple and return it
		return tuple(numericalParts)

	@classmethod
	def checkCMakeVersion(cls):
		"""
		### This function checks if the current installed version, if installed, is above the minimum required version.
		### If no version is provided it just checks if CMake is installed.

		#### Excepts:
		An error message in color red in a string if there is a problem with CMake.
		"""
		# Defining link for cmake installation & update guide
		cmakeInstallURL = 'https://github.com/I2PC/xmipp/wiki/Cmake-update-and-install'

		try:
			# Getting CMake version
			result = subprocess.check_output(["cmake", "--version"]).decode("utf-8")
			cmakVersion = result.split('\n')[0].split()[-1]

			# Checking if installed version is below minimum required
			if CMAKE_MINIMUM_VERSION_V26 and (cls.versionTuple(cmakVersion) < cls.versionTuple(CMAKE_MINIMUM_VERSION_V26)):
				raise Exception(redStr(f"Your CMake version ({cmakVersion}) is below {CMAKE_MINIMUM_VERSION_V26}.\nPlease update your CMake version by following the instructions at {cmakeInstallURL}"))
		except FileNotFoundError:
			raise FileNotFoundError(redStr(f"CMake is not installed.\nPlease install your CMake version by following the instructions at {cmakeInstallURL}"))
		except Exception:
			raise Exception(redStr("Can not get the cmake version.\nPlease visit https://github.com/I2PC/xmipp/wiki/Cmake-troubleshoting"))

	@classmethod
	def parseIndexFile(cls, protocol, indexFile):
		groups, index = {}, 0
		with open(indexFile) as f:
			for line in f:
				if line.startswith('['):
					groups[index] = line.replace('[', '').replace(']', '').strip()
					index += 1
		return groups

	@classmethod
	def translateNamesToIndexGroup(cls, protocol, names):
		"""Translate group name(s) to their numeric index in the index file."""
		indexFile = cls.ensureIndexFile(protocol)
		groups = cls.parseIndexFile(protocol, indexFile)
		invGroups = {v: k for k, v in groups.items()}
		return [invGroups.get(name, name) for name in names]

	@classmethod
	def createIndexFile(cls, protocol, system, inIndex=None, outIndex=None, inputCommands=['q']):
		outIndex = protocol._getExtraPath('indexes.ndx') if not outIndex else outIndex
		outDir = (os.path.dirname(outIndex))
		inIndex = f' -n {os.path.abspath(inIndex)}' if inIndex else ''
		command = f'make_ndx -f {os.path.abspath(system.getSystemFile())}{inIndex} -o {os.path.abspath(outIndex)}'

		if inputCommands[-1] != 'q':
			inputCommands.append('q')
		cls.runGromacsPrintfViewer(printfValues=inputCommands, args=command, cwd=outDir)
		groups = cls.parseIndexFile(protocol, outIndex)
		return outIndex

	@classmethod
	def ensureIndexFile(cls, protocol):
		"""Return the protocol index file, creating it first if it does not exist.
        """
		indexFile = cls.getCustomIndexFile(protocol)
		if os.path.exists(indexFile):
			return os.path.abspath(indexFile)
		indexFile = protocol._getExtraPath('indexes.ndx')
		if os.path.exists(indexFile):
			return os.path.abspath(indexFile)
		inpSystem = protocol.gromacsSystem.get()
		indexFile = inpSystem.getIndexFile()
		if os.path.exists(indexFile):
			return os.path.abspath(indexFile)
		if not os.path.exists(indexFile):
			indexFile = cls.firstIndexCreation(protocol, inpSystem)
		return os.path.abspath(indexFile)

	@classmethod
	def getCustomIndexFile(cls, protocol):
		inputSystem = protocol.gromacsSystem.get()
		inputId = inputSystem.getObjId()
		return protocol.getProject().getTmpPath(f'{inputId}_custom_indexes.ndx')

	@classmethod
	def firstIndexCreation(cls, protocol, groSystem, ligandName=None, modelChains=None, chainLengths=None):
		indexCommands = []

		if ligandName is not None:
			indexCommands.append('1 | 13')
			indexFile = cls.createIndexFile(protocol, groSystem, inputCommands=indexCommands)
		else:
			# Create basic index file with default GROMACS groups
			indexFile = cls.createIndexFile(protocol, groSystem)

		if modelChains is not None and chainLengths is not None and len(modelChains) > 1:
			# Parse existing index file to find the last group number
			groups = cls.parseIndexFile(protocol, indexFile)
			lastGroupIndex = max(groups.keys())

			indexCommands = []
			residuePointer = 1

			for chainId in modelChains:
				chainLength = chainLengths[chainId]
				start = residuePointer
				end = residuePointer + chainLength - 1

				indexCommands.append(f'ri {start}-{end}')

				residuePointer = end + 1

			# Name each newly created group
			for i, chainId in enumerate(modelChains):
				group_number = lastGroupIndex + i + 1
				indexCommands.append(f'name {group_number} chain{chainId}')

			indexCommands.append('q')

			# Create updated index file with per-chain groups
			indexFile = cls.createIndexFile(
				protocol, groSystem,
				inIndex=os.path.abspath(indexFile),
				inputCommands=indexCommands
			)
		return indexFile

	@classmethod
	def checkCudaVersion(cls):
		"""
        ### This function checks if the current installed CUDA version (nvcc) is above the minimum required version.
        ### If no version is provided it just checks if nvcc is installed.
        """
		cudaInstallURL = 'https://developer.nvidia.com/cuda-downloads'
		try:
			# Getting CUDA version from nvcc
			result = subprocess.check_output(["nvcc", "--version"]).decode("utf-8")

			lines = result.split('\n')
			releaseLine = [line for line in lines if 'release' in line][0]
			cudaVersion = releaseLine.split('release ')[1].split(',')[0]

			# Checking if installed version is below minimum required
			if CUDA_MINIMUM_VERSION_V26 and (cls.versionTuple(cudaVersion) < cls.versionTuple(CUDA_MINIMUM_VERSION_V26)):
				raise Exception(redStr(
					f"CUDA version ({cudaVersion}) is below {CUDA_MINIMUM_VERSION_V26}.\nPlease update your CUDA version or set your PATH to a newer version by following instructions at {cudaInstallURL}"))

		except FileNotFoundError:
			raise FileNotFoundError(redStr(
				f"nvcc (CUDA) is not installed or not in your PATH.\n"
				f"Please install CUDA >= {CUDA_MINIMUM_VERSION_V26} or update your PATH by following the instructions at {cudaInstallURL}"))
		except Exception as e:
			raise Exception(redStr(f"Cannot get the CUDA version: {str(e)}\nPlease check your CUDA installation."))

# **************************************************************************
# *
# * Authors:	  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          	  James M. Krieger (jamesmkrieger@gmail.com)
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
import subprocess, multiprocessing
from os.path import join
import os
ENV_YAML_PATH = join(os.path.split(os.path.split(os.path.abspath(__file__))[0])[0],
					 'environment.yaml')

# Scipion em imports
import pwem
from scipion.install.funcs import InstallHelper
from pyworkflow.utils import redStr, yellowStr

# Plugin imports
from .objects import *
from .constants import *

_logo = "gromacs_logo.png"
_references = ['Abraham2015']

class Plugin(pwem.Plugin):
	_name = 'gromacs'
	_version = '1.0'
	_description = 'Gromacs plugin for Scipion-Chem'
	_homeVar = GROMACS_DIC['home']
	_pathVars = [GROMACS_DIC['home'], PLUMED_DIC['home']]
	_supportedVersions = [V2020, V2021, V2024]
	_gromacsName = GROMACS_DIC['name'] + '-' + GROMACS_DIC['version']
	_plumedName = PLUMED_DIC['name'] + '-' + PLUMED_DIC['version']
	_libtorchName = LIBTORCH_DIC['name'] + '-' + LIBTORCH_DIC['version']
	_emmvoxName = EMMIVOX_DIC['name'] + '-' + EMMIVOX_DIC['version']

	@classmethod
	def _defineVariables(cls):
		""" Return and write a variable in the config file. """
		cls._defineEmVar(GROMACS_DIC['home'], GROMACS_DIC['name'] + '-' + GROMACS_DIC['version'])
		cls._defineEmVar(PLUMED_DIC['home'], PLUMED_DIC['name'] + '-' + PLUMED_DIC['version'])
		cls._defineVar(GROMACS_ENV_ACT, f"conda activate {cls.getEnvName()}")

	@classmethod
	def defineBinaries(cls, env):
		""" This function defines all the packages that will be installed. """
		# Checking requirements
		modifiedProcs = cls.checkRequirements(env)

		for ver in LIBTORCH_VERSIONS:
			cls.addLibtorch(env, ver,
				   default=(ver==LIBTORCH_DIC['version']))

		# Installing packages
		for plumed_ver in PLUMED_VERSIONS:
			cls.addPlumed(env, plumed_ver,
				 default=(plumed_ver==PLUMED_DIC['version']))

		for ver in GROMACS_VERSIONS:
			cls.addGromacs(env, modifiedProcs, ver,
				default=(ver==GROMACS_DIC['version']))
		
		for ver in EMMIVOX_VERSIONS:
			cls.addEmmiVox(env, ver,
				  default=(ver==EMMIVOX_DIC['version']))

	@classmethod
	def addLibtorch(cls, env, ver, default=True):
		""" This function downloads Libtorch 2.0.0 for Plumed. """

		# Instantiating install helper
		installer = InstallHelper(LIBTORCH_DIC['name'], join(SCIPION_SOFTWARE, 'em', f"{LIBTORCH_DIC['name']}-{ver}"), ver)

		# Defining some variables
		libTorchFileName = f"{LIBTORCH_DIC['name']}-{LIBTORCH_DIC['version']}.zip"

		ENV_NAME = cls.getEnvName()
		installEnv = [
			cls.getCondaActivationCmd(),
			"if ! conda info --envs | awk '{print $1}' | grep -qx %s; then" % ENV_NAME,
			f'    conda env create -f {ENV_YAML_PATH} -n {ENV_NAME};',
			'fi &&',
			f'conda activate {ENV_NAME} &&',
			'pip install torch==2.0.0+cu118 torchvision==0.15.1+cu118 torchaudio==2.0.1+cu118 ',
			'--index-url https://download.pytorch.org/whl/cu118 &&',
			f'touch ENV_INSTALLED'
		]

		# Installing package
		installer.addCommand(' '.join(installEnv), 'ENV_INSTALLED')\
			.getExtraFile(cls._getLibTorchDownloadUrl(LIBTORCH_DIC['version']), 'LIBTORCH_DOWNLOADED', fileName=libTorchFileName)\
			.addCommand(f'unzip {libTorchFileName}', 'LIBTORCH_EXTRACTED')\
			.addPackage(env, dependencies=['wget', 'unzip'], default=default)

	@classmethod
	def addPlumed(cls, env, ver, default=True):
		""" This function installs Plumed's package. """

		# Instantiating install helper
		plumedLocation = join(SCIPION_SOFTWARE, 'em', f"{PLUMED_DIC['name']}-{ver}")
		libtorchLocation = cls._getLocation(LIBTORCH_DIC, marker='LIBTORCH_EXTRACTED')
		installer = InstallHelper(PLUMED_DIC['name'], plumedLocation, ver)

		# Defining some variables
		plumedFileName = f"{PLUMED_DIC['name']}-{ver}.tar.gz"

		ENV_NAME = cls.getEnvName()
		installPlumed = [
			cls.getCondaActivationCmd(),
			f'conda activate {ENV_NAME} &&'
		]

		if cls._isInstalled(LIBTORCH_DIC, marker='LIBTORCH_EXTRACTED', location=libtorchLocation):
			libtorchLocation += '/libtorch'
			installPlumed.extend([
				f'CPPFLAGS="-I{libtorchLocation}/include -I{libtorchLocation}/include/torch/csrc/api/include"',
				f'LDFLAGS="-L{libtorchLocation}/lib -Wl,-rpath,{libtorchLocation}/lib"',
				'./configure --enable-libtorch --enable-modules=pytorch',
				f"--prefix={plumedLocation}/install &&",
				'make && make install &&',
				f'touch PLUMED_INSTALLED'
			])
		else:
			installPlumed.extend([
				f"./configure --prefix={plumedLocation}/install &&",
				'make && make install &&',
				f'touch PLUMED_INSTALLED'
			])

		# Installing package
		installer.getExtraFile(cls._getPlumedDownloadUrl(ver), 'PLUMED_DOWNLOADED', fileName=plumedFileName)\
			.addCommand(f'tar -xf {plumedFileName} --strip-components 1', 'PLUMED_EXTRACTED')\
			.addCommand(' '.join(installPlumed), 'PLUMED_INSTALLED')\
			.addPackage(env, dependencies=['wget', 'tar', 'cmake', 'make'], default=default)

	@classmethod
	def addGromacs(cls, env, modifiedProcs, ver, default=True):
		""" This function installs Gromacs's package. """
		# Target suffix for MPI installation
		mpiExt = '_MPI'
		plumedExt = '_PLUMED'

		# Instantiating install helper
		installer = InstallHelper(GROMACS_DIC['name'], join(SCIPION_SOFTWARE, 'em', f"{GROMACS_DIC['name']}-{ver}"), ver)

		# Defining some variables
		gromacsFileName = f"gromacs-{ver}.tar.gz"
		
		charmInnerLocation = join('share', 'top')
		charmFileName = 'charmm36-feb2021.ff.tgz'

		plumed_ver = cls._getInstalledVersion(PLUMED_DIC)
		plumedLocation = join(SCIPION_SOFTWARE, 'em', f"{PLUMED_DIC['name']}-{plumed_ver}")

		normalInnerLocation = 'build'
		mpiInnerLocation = 'build_mpi'
		plumedMpiInnerLocation = 'build_mpi_plumed'

		# If number of processors has been modified, show message
		if modifiedProcs:
			message = "\\nWARNING: Only 1 process has been defined to install Gromacs.\\n"
			message += f"This will take a very long time. Instead, the number of parallel processes has been changed to the maximum available in your system: {env.getProcessors()}."
			installer.addCommand(f'echo -e "{yellowStr(message)}"', 'WARNING_SHOWN')

		ENV_NAME = cls.getEnvName()

		patchGromacsWithPlumed = [
			cls.getCondaActivationCmd(),
			f'conda activate {ENV_NAME} &&',
			f"export PATH=$PATH:{plumedLocation}/bin &&",
			f"export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{plumedLocation}/lib &&",
			f"export PLUMED_KERNEL={plumedLocation}/lib/libplumedKernel.so &&",			
			f'echo "{PATCH_DIC.get(plumed_ver, None).get(ver, None)}" >> patch_option.txt',
			'plumed patch -p --runtime < patch_option.txt'
		]

		# Installing package
		installer.getExtraFile(cls._getGromacsDownloadUrl(ver), 'GROMACS_DOWNLOADED', fileName=gromacsFileName)\
			.addCommand(f'tar -xf {gromacsFileName} --strip-components 1', 'GROMACS_EXTRACTED')

		CUDA_ARCH_FLAG = '-DGMX_CUDA_TARGET_SM="50;52;60;61;70;75;80"' if ver==V2021 else '-DCUDA_ARCH_BIN=all'

		installer.getExtraFile('http://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/charmm36-feb2021.ff.tgz', 'CHARM_DOWNLOADED', location=charmInnerLocation, fileName=charmFileName)\
			.addCommand(f'tar -xf {charmFileName}', 'CHARM_EXTRACTED', workDir=charmInnerLocation)\
			.addCommand(f'mkdir {normalInnerLocation} {mpiInnerLocation}', 'BUILD_DIRS_MADE')\
			.addCommand(f'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA {CUDA_ARCH_FLAG} -DCMAKE_INSTALL_PREFIX={cls._getLocation(GROMACS_DIC, ver)}/install -DGMX_FFT_LIBRARY=fftw3',
		   				'GROMACS_BUILT', workDir=normalInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()}', 'GROMACS_COMPILED', workDir=normalInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()} install', 'GROMACS_INSTALLED', workDir=normalInnerLocation)\
			.addCommand(f'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA {CUDA_ARCH_FLAG} -DCMAKE_INSTALL_PREFIX={cls._getLocation(GROMACS_DIC, ver)}/install{mpiExt.lower()} -DGMX_FFT_LIBRARY=fftw3 -DGMX_MPI=ON',
		   				'GROMACS_BUILT' + mpiExt, workDir=mpiInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()}', 'GROMACS_COMPILED' + mpiExt, workDir=mpiInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()} install', 'GROMACS_INSTALLED' + mpiExt, workDir=mpiInnerLocation)

		if (cls._isInstalled(PLUMED_DIC, marker='PLUMED_INSTALLED', location=plumedLocation) and
			PATCH_DIC.get(plumed_ver, None).get(ver, None) is not None):
			installer.addCommand(' '.join(patchGromacsWithPlumed), 'PLUMED_PATCHED')\
			.addCommand(f'mkdir {plumedMpiInnerLocation}', 'PLUMED_BUILD_DIR_MADE')\
			.addCommand(f'cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA {CUDA_ARCH_FLAG} -DCMAKE_INSTALL_PREFIX={cls._getLocation(GROMACS_DIC, ver)}/install{mpiExt.lower()} -DGMX_FFT_LIBRARY=fftw3 -DGMX_MPI=ON',
						'GROMACS_BUILT' + mpiExt + plumedExt, workDir=plumedMpiInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()}', 'GROMACS_COMPILED' + mpiExt + plumedExt, workDir=plumedMpiInnerLocation)\
			.addCommand(f'make -j{env.getProcessors()} install', 'GROMACS_INSTALLED' + mpiExt + plumedExt, workDir=plumedMpiInnerLocation)

		installer.addPackage(env, dependencies=['wget', 'tar', 'cmake', 'make'], default=default)

	@classmethod
	def addEmmiVox(cls, env, ver, default=True):
		""" This function downloads EMMIVox 0.1 scripts for use with Gromacs and Plumed. """

		# Instantiating install helper
		installer = InstallHelper(EMMIVOX_DIC['name'], join(SCIPION_SOFTWARE, 'em', f"{EMMIVOX_DIC['name']}-{ver}"), ver)

		# Defining some variables
		emmiVoxFileName = f"{EMMIVOX_DIC['name']}-{EMMIVOX_DIC['version']}"
		mapFileName1 = 'emd_13223_half_map_1.map.gz'
		mapFileName2 = 'emd_13223_half_map_2.map.gz'
		mapFileNameMain = 'emd_13223.map.gz'

		# Installing package
		installer.getCloneCommand(cls._getEmmiVoxDownloadUrl(), targeName='EMMIVOX_CLONED', binaryFolderName=emmiVoxFileName)\
			.addCommand(f'cd {emmiVoxFileName} && git checkout scipion && cd ..', 'EMMIVOX_SCIPION_CHECKED_OUT')\
			.getExtraFile(cls._getEmmiTutorialHalfMapUrl(1), 'EMMIVOX_HALF_MAP_1_DOWNLOADED', fileName=mapFileName1)\
			.addCommand(f'gunzip {mapFileName1}', 'EMMIVOX_HALF_MAP_1_EXTRACTED')\
			.getExtraFile(cls._getEmmiTutorialHalfMapUrl(2), 'EMMIVOX_HALF_MAP_2_DOWNLOADED', fileName=mapFileName2)\
			.addCommand(f'gunzip {mapFileName2}', 'EMMIVOX_HALF_MAP_2_EXTRACTED')\
			.getExtraFile(cls._getEmmiTutorialMainMapUrl(), 'EMMIVOX_MAIN_MAP_DOWNLOADED', fileName=mapFileNameMain)\
			.addCommand(f'gunzip {mapFileNameMain}', 'EMMIVOX_MAIN_MAP_EXTRACTED')\
			.addPackage(env, dependencies=['wget', 'unzip'], default=default)

	@classmethod
	def runGromacs(cls, protocol, program='gmx', args='', cwd=None, mpi=False, **kwargs):
		""" Run Gromacs command from a given protocol. """
		protocol.runJob(cls.getGromacsBin(program, mpi=mpi), args, cwd=cwd, **kwargs)

	@classmethod
	def runGromacsPrintf(cls, printfValues, args, cwd, mpi=False):
		""" Run Gromacs command from a given protocol. """
		printfValues = list(map(str, printfValues))
		program = 'printf "{}\n" | {} '.format('\n'.join(printfValues), cls.getGromacsBin(mpi=mpi))
		print('Running: ', program, args)
		subprocess.check_call(program + args, cwd=cwd, shell=True)

	@classmethod
	def getGromacsBin(cls, program='gmx', mpi=False):
		mpiExt = '_mpi' if mpi else ''
		return join(cls._getLocation(GROMACS_DIC, marker=f'GROMACS_INSTALLED{mpiExt.upper()}'),
			  f'install{mpiExt}/bin/{program}{mpiExt}')

	@classmethod
	def getSourceGromacsCmd(cls, mpi=False):
		mpiExt = '_mpi' if mpi else ''
		return ". " + join(
			cls._getLocation(GROMACS_DIC, marker=f'GROMACS_INSTALLED{mpiExt.upper()}'),
			f'install{mpiExt}/bin/GMXRC')

	@classmethod  # Test that
	def getEnviron(cls):
		pass

	@classmethod
	def runPlumed(cls, protocol, program='plumed', args='', cwd=None, **kwargs):
		""" Run Plumed command from a given protocol. """
		protocol.runJob(cls.getPlumedBin(program), args, cwd=cwd, **kwargs)

	@classmethod
	def runPlumedPrintf(cls, printfValues, args, cwd, program='plumed'):
		""" Run Plumed command from a given protocol. """
		printfValues = list(map(str, printfValues))
		program = 'printf "{}\n" | {} '.format('\n'.join(printfValues),
										 cls.getPlumedBin(program))
		print('Running: ', program, args)
		subprocess.check_call(program + args, cwd=cwd, shell=True)

	@classmethod
	def getPlumedBin(cls, program='plumed', location=None):
		if location is None:
			plumed_ver = cls._getInstalledVersion(PLUMED_DIC)
			location = cls._getLocation(PLUMED_DIC, plumed_ver)
		return f"{location}/install/bin/{program}"

	@classmethod
	def _getLocation(cls, dic, version=None, marker=None):
		if version is None:
			version = cls._getInstalledVersion(dic, marker=marker)
		return join(SCIPION_SOFTWARE, 'em', f"{dic['name']}-{version}")
	
	@classmethod
	def _getLibTorchLocation(cls):
		return cls._getLocation(LIBTORCH_DIC, marker="LIBTORCH_EXTRACTED")
	
	@classmethod
	def _getEmmiVoxLocation(cls):
		return join(cls._getLocation(EMMIVOX_DIC, marker="EMMIVOX_CLONED"),
			  EMMIVOX_DIC['name'] + "-" + cls._getInstalledVersion(
				  EMMIVOX_DIC, marker='EMMIVOX_CLONED'))

	@classmethod
	def _getEmmiVoxTutorialLocation(cls):
		return join(cls._getEmmiVoxLocation(), "tutorials")

	@classmethod
	def _getEmmiVoxScriptLocation(cls):
		return join(cls._getEmmiVoxLocation(), "scripts")

	@classmethod
	def getEmmiVoxProgram(cls, program='', python=False, location=None):
		""" Create EMMIVox script command line. 
		
		:arg python: whether it needs to be run with Python
			Otherwise, use bash. Default is **False**
		:type python: bool
		"""
		baseLocation = cls._getEmmiVoxLocation()
		if location is None:
			# assume it's in the scripts directory
			location = join(baseLocation, "scripts")

		elif not os.path.exists(location) and os.path.exists(join(baseLocation, location)):
			location = join(baseLocation, location)

		if python is False and program.endswith('.py'):
			python=True

		if python:
			fullProgram = '%s %s && python %s' % (
				cls.getCondaActivationCmd(), cls.getEnvName(),
				join(location, program))
		else:
			fullProgram = '%s && bash %s' % (
				cls.getSourceGromacsCmd(mpi=True),
				join(location, program))

		return fullProgram

	@classmethod
	def _getLibTorchDownloadUrl(cls, ver):
		return 'https://download.pytorch.org/libtorch/cu118/libtorch-cxx11-abi-shared-with-deps-{}%2Bcu118.zip'.format(ver)

	@classmethod
	def _getPlumedDownloadUrl(cls, ver):
		if ver == MASTER:
			return 'https://github.com/plumed/plumed2/archive/refs/heads/{}.tar.gz'.format(ver)
		return 'https://github.com/plumed/plumed2/archive/refs/tags/v{}.tar.gz'.format(ver)

	@classmethod
	def _getGromacsDownloadUrl(cls, ver):
		return 'https://ftp.gromacs.org/gromacs/gromacs-{}.tar.gz'.format(ver)
	
	@classmethod
	def _getEmmiVoxDownloadUrl(cls):
		return 'git@github.com:jamesmkrieger/EMMIVox.git'
	
	@classmethod
	def _getEmmiTutorialHalfMapUrl(cls, mapNum=1):
		return f'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-13223/other/emd_13223_half_map_{mapNum}.map.gz'
	
	@classmethod
	def _getEmmiTutorialMainMapUrl(cls):
		return 'https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-13223/map/emd_13223.map.gz'

	# ---------------------------------- Utils functions  -----------------------
	@classmethod
	def checkRequirements(cls, env):
		""" This function checks if the software requirements are being met. """
		cls.checkCMakeVersion()
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
		cmakeInstallURL = 'https://github.com/I2PC/xmipp/wiki/Cmake-update-and/install'

		try:
			# Getting CMake version
			result = subprocess.check_output(["cmake", "--version"]).decode("utf-8")
			cmakVersion = result.split('\n')[0].split()[-1]

			# Checking if installed version is below minimum required
			if CMAKE_MINIMUM_VERSION and (cls.versionTuple(cmakVersion) < cls.versionTuple(CMAKE_MINIMUM_VERSION)):
				raise Exception(redStr(f"Your CMake version ({cmakVersion}) is below {CMAKE_MINIMUM_VERSION}.\nPlease update your CMake version by following the instructions at {cmakeInstallURL}"))
		except FileNotFoundError:
			raise FileNotFoundError(redStr(f"CMake is not installed.\nPlease install your CMake version by following the instructions at {cmakeInstallURL}"))
		except Exception:
			raise Exception(redStr("Can not get the cmake version.\nPlease visit https://github.com/I2PC/xmipp/wiki/Cmake-troubleshoting"))

	@classmethod
	def _getNvccCudaVersion(cls):
		"""
		### This function checks what version of cuda and nvcc is active and returns the result

		### Excepts:
		An error message in color red in a string if there is a problem with nvcc
		"""
		try:
			# Getting CMake version
			result = subprocess.check_output(["nvcc", "--version"]).decode("utf-8")
			nvccVersion = result.split('\n')[-2].split('_')[1].split('/')[0][:4]
			return nvccVersion

			# Checking if installed version is below minimum required
		except FileNotFoundError:
			raise FileNotFoundError(redStr(f"nvcc is not installed.\nPlease install A CUDA toolkit in development to include nvcc."))
		except Exception:
			raise Exception(redStr("Can not get the nvcc version.\nPlease install A CUDA toolkit in development to include nvcc."))

	@classmethod
	def _isInstalled(cls, dic, marker=None, location=None):
		"""Check if package is installed by looking for its marker file."""
		if location is None:
			home = cls.getVar(dic['home'])
			if home is None:
				return False
		else:
			home = location

		if marker is None:
			marker = dic['name'].upper()+'_INSTALLED'
		markerFile = os.path.join(home, marker)
		return os.path.exists(markerFile)

	@classmethod
	def _getInstalledVersion(cls, dic, marker=None):
		"""Return installed version from directory name or version file."""
		active = cls._getActiveVersion(dic)
		all_versions = cls._getInstalledVersions(dic, marker=marker)
		if active is not None:
			return active
		elif len(all_versions):
			return all_versions[-1]
		return None

	@classmethod
	def _getInstalledVersions(cls, dic, marker=None):
		"""
		Return a list of installed versions for a software dictionary,
		using _isInstalled() to verify marker files.
		"""
		em_dir = os.path.join(SCIPION_SOFTWARE, "em")
		installed_versions = []

		if not os.path.exists(em_dir):
			return installed_versions

		for entry in os.listdir(em_dir):
			path = os.path.join(em_dir, entry)
			if os.path.isdir(path):
				if cls._isInstalled(dic, marker=marker, location=path):
					# parse version from directory name (last part after '-')
					version = entry.split('-')[-1]
					installed_versions.append(version)

		return installed_versions

	@classmethod
	def _getActiveVersion(cls, dic):
		"""
		Return the version currently referenced by dic['home'],
		or None if home is not set or doesn't exist.
		"""
		home = cls.getVar(dic.get('home', None))
		if home and os.path.exists(home):
			return os.path.basename(home).split('-')[-1]
		return None

	@classmethod
	def getEnvName(cls):
		return LIBTORCH_DIC['name']+"-"+LIBTORCH_DIC['version']

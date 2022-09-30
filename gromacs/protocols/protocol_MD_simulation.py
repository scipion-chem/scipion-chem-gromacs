# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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


"""
This module will perform energy minimizations for the system
"""
import os, glob, shutil

from pyworkflow.protocol import params
from pyworkflow.utils import Message, runJob, createLink
from pwem.protocols import EMProtocol

from pwchem.utils import natural_sort

from gromacs.objects import *
from gromacs.constants import *
from gromacs import Plugin as gromacsPlugin

from multiprocessing import cpu_count

class GromacsMDSimulation(EMProtocol):
    """
    This protocol will perform energy minimization on the system previosly prepared by the protocol "system prepartion".
    This step is necessary to energy minize the system in order to avoid unwanted conformations.
    """
    _label = 'Molecular dynamics simulation'
    _ensemTypes = ['Energy min', 'NVT',  'NPT']

    _integrators = ['steep', 'cg']
    _thermostats = ['no', 'Berendsen', 'Nose-Hoover', 'Andersen', 'Andersen-massive', 'V-rescale']
    _barostats = ['no', 'Berendsen', 'Parrinello-Rahman']
    #_coupleStyle = ['isotropic', 'semiisotropic', 'anisotropic'] #check
    _restraintTypes = ['None', 'Protein', 'Protein-H', 'MainChain', 'BackBone', 'C-alpha']

    _paramNames = ['simTime', 'timeStep', 'nStepsMin', 'emStep', 'emTol', 'timeNeigh', 'saveTrj', 'trajInterval',
                   'temperature', 'tempRelaxCons', 'tempCouple', 'pressure', 'presRelaxCons', 'presCouple',
                   'restraintForce']
    _enumParamNames = ['integrator', 'ensemType', 'thermostat', 'barostat', 'restraints']
    _defParams = {'simTime': 100, 'timeStep': 0.002, 'nStepsMin': 50000, 'emStep': 0.002, 'emTol': 1000.0,
                  'timeNeigh': 10, 'saveTrj': False, 'trajInterval': 1.0,
                  'temperature': 300.0, 'tempRelaxCons': 0.1, 'tempCouple': -1, 'integrator': 'cg',
                  'pressure': 1.0, 'presRelaxCons': 2.0, 'presCouple': -1,'restraintForce': 50.0,
                  'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman',
                  'restraints': 'None'}

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
      EMProtocol.__init__(self, **kwargs)


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=0)

        form.addSection(label=Message.LABEL_INPUT)
        form.addHidden(params.USE_GPU, params.BooleanParam,
                      label='Use GPU for execution',
                      default=False,
                      help="GPU may have several cores. Set it one if "
                           "you don't know what we are talking about but you have a GPU."
                           "For DARC, first core index is 1, second 2, and so on. Write 0 if you do not want"
                           "to use GPU")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=params.LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used (Comma separated)")

        form.addParam('gromacsSystem', params.PointerParam, label="Input Gromacs System: ",
                      pointerClass='GromacsSystem',
                      help='Gromacs solvated system to be simulated')
        form.addParam('prevTrj', params.BooleanParam, default=False,
                      label="Concatenate with previous trajectory: ", expertLevel=params.LEVEL_ADVANCED,
                      help='Include the trajectory stored in the input system (if found). \n'
                           'It will only be included if all the stages in this protocol have their trajectories saved,'
                           ' for the sake of continuity')
        form.addParam('cptTime', params.FloatParam, default=15,
                      expertLevel=params.LEVEL_ADVANCED,
                      label='Checkpoint time (min):',
                      help='Time of execution for saving a checkpoint. Checkpoints allow the simulation to be continued'
                           ' from that moment, instead of restarting from the start. \nIn scipion, the simulation will'
                           'be continued from the checkpoint if the option "Continue" is chosen after the protocol was'
                           'stopped for any reason')

        group = form.addGroup('Ensemble')
        group.addParam('ensemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
                       help='Type of simulation to perform in the step: Energy minimization, NVT or NPT\n'
                            'https://manual.gromacs.org/5.1.1/user-guide/mdp-options.html')

        group.addParam('integrator', params.EnumParam,
                      label='Simulation integrator: ', condition='ensemType==0',
                      choices=self._integrators, default=0,
                      help='Type of integrator to use in simulation.')

        line = group.addLine('Temperature settings: ', condition='ensemType!=0',
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Relaxation time constant for thermostat (ps)')
        line.addParam('temperature', params.FloatParam, default=300, condition='ensemType!=0',
                      label='Temperature: ')
        line.addParam('thermostat', params.EnumParam, default=5, condition='ensemType!=0',
                      label='Thermostat: ', choices=self._thermostats)
        line.addParam('tempRelaxCons', params.FloatParam, default=0.1,
                      label='Temperature constant (ps)[tau-t]: ', expertLevel=params.LEVEL_ADVANCED)
        line.addParam('tempCouple', params.IntParam, default=-1,
                      label='Coupling frequency [nsttcouple]: ', expertLevel=params.LEVEL_ADVANCED)

        line = group.addLine('Pressure settings: ', condition='ensemType==2',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('pressure', params.FloatParam, default=1.0,
                      label='   Pressure (bar):   ')
        line.addParam('barostat', params.EnumParam, default=2,
                      label='  Barostat type:   ', choices=self._barostats)
        line.addParam('presRelaxCons', params.FloatParam, default=2.0,
                      label='   Pressure constant (ps)[tau-p]:   ', expertLevel=params.LEVEL_ADVANCED)
        line.addParam('presCouple', params.IntParam, default=-1,
                      label='Coupling frequency [nstpcouple]: ', expertLevel=params.LEVEL_ADVANCED)
        #group.addParam('coupleStyle', params.EnumParam, default=0, condition='ensemType==2',
        #               label='Pressure coupling style: ', choices=self._coupleStyle,
        #               expertLevel=params.LEVEL_ADVANCED)

        group = form.addGroup('Trajectory', condition='ensemType!=0')
        group.addParam('saveTrj', params.BooleanParam, default=self._defParams['saveTrj'],
                       label="Save trajectory: ", condition='ensemType!=0',
                       help='Save trajectory of the atoms during stage simulation.'
                            'The output will concatenate those trajectories which appear after the last stage '
                            'where the trajectory was not saved.')
        group.addParam('trajInterval', params.FloatParam, default=self._defParams['trajInterval'],
                       label='Interval time (ps):', condition='ensemType!=0 and saveTrj',
                       help='Time between each frame recorded in the simulation (ps)')

        group = form.addGroup('Simulation time')
        group.addParam('simTime', params.FloatParam, default=100,
                       label='Simulation time (ps):', condition='ensemType!=0',
                       help='Total time of the simulation stage (ps)')
        group.addParam('timeStep', params.FloatParam, default=0.002,
                       label='Simulation time steps (ps)[dt]:', condition='ensemType!=0',
                       help='Time of the steps for simulation (ps)[dt]')

        group.addParam('nStepsMin', params.IntParam, default=50000,
                       label='Maximum number of steps:', condition='ensemType==0',
                       help='Maximum number of (minimization) steps to perform')
        group.addParam('emTol', params.FloatParam, default=1000.0,
                       label='Maximum force objective:', condition='ensemType==0',
                       help='Stop minimization when the maximum force < x kJ/mol/nm.\n'
                            'https://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html#mdp-emtol')
        group.addParam('emStep', params.FloatParam, default=0.002,
                       label='Initial step-size (nm)[emstep]:', condition='ensemType==0',
                       help='Initial step-size (nm)[emstep].\n'
                            'https://manual.gromacs.org/documentation/2018/user-guide/mdp-options.html#mdp-emstep')

        group.addParam('timeNeigh', params.IntParam, default=10,
                       label='Frequency to update the neighbor list (steps)[nstlist]:',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Frequency to update the neighbor list (and the long-range forces, when using twin-range '
                            'cut-offs). When this is 0, the neighbor list is made only once. With energy minimization '
                            'the neighborlist will be updated for every energy evaluation when nstlist is '
                            'greater than 0.')

        group = form.addGroup('Restraints')
        group.addParam('restraints', params.EnumParam, default=0,
                       label='Restraints: ', choices=self._restraintTypes,
                       help='Restraint movement of specific groups of atoms')
        group.addParam('restraintForce', params.FloatParam, default=50,
                       label='Restraint force constant: ', condition='restraints!=0',
                       help='Restraint force applied to the selection (kcal/mol/Å2)')

        group = form.addGroup('Summary')
        group.addParam('insertStep', params.StringParam, default='',
                       label='Insert relaxation step number: ',
                       help='Insert the defined relaxation step into the workflow on the defined position.\n'
                            'The default (when empty) is the last position')
        group.addParam('summarySteps', params.TextParam, width=120, readOnly=True,
                       label='Summary of steps',
                       help='Summary of the defined steps. \nManual modification will have no '
                            'effect, use the wizards to add / delete the steps')
        group.addParam('deleteStep', params.StringParam, default='',
                       label='Delete relaxation step number: ',
                       help='Delete the step of the specified index from the workflow.')
        group.addParam('watchStep', params.StringParam, default='',
                       label='Watch relaxation step number: ',
                       help='''Watch the parameters step of the specified index from the workflow..\n
                               This might be useful if you want to change some parameters of a predefined step.\n
                               However, the parameters are not changed until you add the new step (and probably\n
                               you may want to delete the previous unchanged step)''')
        group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self.createGUISummary()
        i = 1
        for wStep in self.workFlowSteps.get().strip().split('\n'):
            self._insertFunctionStep('simulateStageStep', wStep, i)
            i += 1
        self._insertFunctionStep('createOutputStep')

    def simulateStageStep(self, wStep, i):
      if wStep in ['', None]:
          msjDic = self.createMSJDic()
      else:
          msjDic = eval(wStep)

      mdpFile = self.generateMDPFile(msjDic, str(i))
      tprFile = self.callGROMPP(mdpFile)
      saveTrj = msjDic['saveTrj'] if msjDic['ensemType'] != 'Energy min' else False

      self.callMDRun(tprFile, saveTrj=saveTrj)

    def createOutputStep(self):
        lastGroFile, lastTopoFile, lastTprFile = self.getPrevFinishedStageFiles()
        if self.gromacsSystem.get().hasTrajectory() and self.prevTrj.get():
            oriGroFile = self.gromacsSystem.get().getOriStructFile()
        else:
            oriGroFile = self.gromacsSystem.get().getSystemFile()

        localGroFile, localTopFile = self._getPath('outputSystem.gro'), self._getPath('systemTopology.top')
        shutil.copyfile(lastGroFile, localGroFile)
        shutil.copyfile(lastTopoFile, localTopFile)

        outTrj = self.concatTrjFiles(outTrj='outputTrajectory.xtc', tprFile=lastTprFile)

        outSystem = GromacsSystem(filename=localGroFile, oriStructFile=oriGroFile,
                                  tprFile=lastTprFile)
        outSystem.setTopologyFile(localTopFile)
        outSystem.setChainNames(self.gromacsSystem.get().getChainNames())
        if outTrj:
            outSystem.setTrajectoryFile(outTrj)
            outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())

        self._defineOutputs(outputSystem=outSystem)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
      fnSummary = self._getExtraPath("summary.txt")
      if not os.path.exists(fnSummary):
        summary = ["No summary information yet."]
      else:
        fhSummary = open(fnSummary, "r")
        summary = []
        for line in fhSummary.readlines():
          summary.append(line.rstrip())
        fhSummary.close()
      return summary

    def writeSummaryLine(self, msjDic):
        ensemType = msjDic['ensemType']
        if ensemType == 'Energy min':
            sumStr = 'Minimization ({}): {} steps, {} objective force'.format(msjDic['integrator'], msjDic['nStepsMin'],
                                                                              msjDic['emTol'])
        else:
            sumStr = 'MD simulation: {} ps, {} ensemble'.format(msjDic['simTime'], ensemType)
  
        if msjDic['restraints'] != 'None':
          sumStr += ', restraint on {}'.format(msjDic['restraints'])
        sumStr += ', {} K\n'.format(msjDic['temperature'])
        return sumStr

    def createSummary(self, msjDic=None):
        '''Creates the displayed summary from the internal state of the steps'''
        if not msjDic:
            sumStr = ''
            for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
                if dicLine != '':
                    msjDic = eval(dicLine)
                    msjDic = self.addDefaultForMissing(msjDic)
                    sumStr += '{}) {}'.format(i+1, self.writeSummaryLine(msjDic))
        else:
            msjDic = self.addDefaultForMissing(msjDic)
            sumStr = self.writeSummaryLine(msjDic)
        return sumStr

    def createGUISummary(self):
        with open(self._getExtraPath("summary.txt"), 'w') as f:
            if self.workFlowSteps.get():
                f.write(self.createSummary())
            else:
                f.write(self.createSummary(self.createMSJDic()))

    def _methods(self):
        methods = []

        if self.isFinished():
            methods.append('The methods used to perform the energy minimization protocol have been *"gmx grompp"* to '
                           'create the tpr files and *"gmx mdrun"* to run the minimization simulation.')

        return methods

    def _warnings(self):
        warns = []
        #Global warnings
        inpSystem = self.gromacsSystem.get()
        if str(inpSystem.getForceField()).startswith('gromos'):
            warns.append('\nStep all : GROMOS force field is not recommended by GROMACS: '
                         'https://chemrxiv.org/engage/chemrxiv/article-details/60c74701bdbb895afaa38ce2')
        elif str(inpSystem.getForceField()).startswith('charmm36') and 'CA' in inpSystem.getIons() and \
                inpSystem.getIons()['CA'] > 20:
          warns.append('\nStep all : More than 20 non-matching atom names because of (Cal - CA) in charmm36. '
                       'Cal from topology file will be used')

        prevTrj = False
        for step, wStep in enumerate(self.workFlowSteps.get().strip().split('\n')):
            if wStep in ['', None]:
                msjDic = self.createMSJDic()
            else:
                msjDic = eval(wStep)

            if msjDic['ensemType'] == 'NPT' and msjDic['barostat'] != 'Berendsen' and not prevTrj and msjDic['saveTrj']:
                warns.append('\nStep {} : Berendsen is the barostat recommended for system equilibration, '
                             '{} might not be the best option for the first trajectory saved'.
                             format(step+1, msjDic['barostat']))
            if msjDic['saveTrj']:
                prevTrj = True

            if msjDic['ensemType'] != 'Energy min':
                tCoup = msjDic['timeNeigh'] if msjDic['tempCouple'] == -1 else msjDic['tempCouple']
                if msjDic['thermostat'] == 'Nose-Hoover' and 20*tCoup*msjDic['timeStep'] > msjDic['tempRelaxCons']:
                  warns.append('\nStep {} : For proper integration of the Nose-Hoover thermostat, tau-t ({}) should '
                               'be at least 20 times larger than nsttcouple*dt ({}*{})'.format(step+1,
                                msjDic['tempRelaxCons'], tCoup, msjDic['timeStep']))

        return warns

    def _validate(self):
        vals = []
        for step, wStep in enumerate(self.workFlowSteps.get().strip().split('\n')):
            if wStep in ['', None]:
                msjDic = self.createMSJDic()
            else:
                msjDic = eval(wStep)
            if msjDic['ensemType'] != 'Energy min':
              if 'Andersen' in msjDic['thermostat'] and msjDic['integrator'] == 'md':
                  vals.append('Step {} : Andersen temperature control not supported for integrator md.'.format(step+1))
        return vals

######################## UTILS ##################################

    def countSteps(self):
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

    def createMSJDic(self):
        msjDic = {}
        for pName in self._paramNames:
            if hasattr(self, pName):
                msjDic[pName] = getattr(self, pName).get()
            else:
                print('Something is wrong with parameter ', pName)

        for pName in self._enumParamNames:
            if hasattr(self, pName):
                msjDic[pName] = self.getEnumText(pName)
            else:
                print('Something is wrong with parameter ', pName)
        return msjDic

    def addDefaultForMissing(self, msjDic):
        '''Add default values for missing parameters in the msjDic'''
        for pName in [*self._paramNames, *self._enumParamNames]:
            if not pName in msjDic:
                msjDic[pName] = self._defParams[pName]
        return msjDic

    def generateMDPFile(self, msjDic, mdpStage):
        stageDir = self._getExtraPath('stage_{}'.format(mdpStage))
        if not os.path.exists(stageDir):
            os.mkdir(stageDir)

        mdpFile = os.path.join(stageDir, 'stage_{}.mdp'.format(mdpStage))
        if os.path.exists(mdpFile): return mdpFile

        restr = msjDic['restraints']
        if restr != 'None':
            rSuffix = msjDic['restraints'] + '_stg%s' % mdpStage
            self.gromacsSystem.get().defineNewRestriction(index=msjDic['restraints'], energy=msjDic['restraintForce'],
                                                          restraintSuffix=rSuffix, outDir=stageDir)
            restrStr = RESTR_STR.format(rSuffix.upper())
        else:
            restrStr = ''

        neighlist = msjDic['timeNeigh']
        if msjDic['ensemType'] == 'Energy min':
            emStep, nSteps, emTol = msjDic['emStep'], msjDic['nStepsMin'], msjDic['emTol']
            integStr = msjDic['integrator']
            tStepsStr = TSTEP_EM.format(emTol, emStep)
            dispCorrStr, outControlStr = '', ''
            bondParStr, electroStr = '', ''
            tempStr, presStr = '', ''
            velStr, velParStr = 'no', ''
            nTrj = round(msjDic['trajInterval'] / emStep)
        else:
            tSteps = msjDic['timeStep']
            nSteps = round(msjDic['simTime'] / tSteps)
            integStr = 'md'
            tStepsStr = TSTEP_EQ.format(tSteps)
            dispCorrStr = DISP_CORR
            electroStr = ELECTROSTATICS
            nstcomm = 10 if msjDic['thermostat'] != 'Andersen' else 1
            tempStr = TEMP_SETTING.format(msjDic['thermostat'],
                                          *[msjDic['temperature']]*2, *[msjDic['tempRelaxCons']]*2,
                                          msjDic['tempCouple'], nstcomm)

            if msjDic['ensemType'] == 'NVT':
                presStr = ''
            elif msjDic['ensemType'] == 'NPT':
                presStr = PRES_SETTING.format(msjDic['barostat'], 'isotropic',
                                              msjDic['pressure'], msjDic['presRelaxCons'],
                                              msjDic['presCouple'])

            if self.checkIfPrevTrj(mdpStage):
                bondParStr = BONDED_PARAMS.format('yes')
                velStr, velParStr = 'no', ''
            else:
                bondParStr = BONDED_PARAMS.format('no')
                velStr = 'yes'
                velParStr = VEL_GEN.format(msjDic['temperature'])
            nTrj = round(msjDic['trajInterval'] / tSteps)

        if msjDic['saveTrj']:
            outControlStr = OUTPUT_CONTROL.format(*[nTrj]*3)
            print('Number of frames: ', nSteps/nTrj)
        else:
            outControlStr = OUTPUT_CONTROL.format(*[0]*2, nTrj)

        title = 'Stage {}: {}, {} ps'.format(mdpStage, msjDic['ensemType'], msjDic['simTime'])
        mdpStr = MDP_STR.format(restrStr, integStr, nSteps, tStepsStr, neighlist, dispCorrStr, outControlStr,
                                bondParStr, electroStr, tempStr, presStr, velStr, velParStr)

        with open(mdpFile, 'w') as f:
            f.write(mdpStr)

        return mdpFile

    def callGROMPP(self, mdpFile):
        stageDir = os.path.dirname(mdpFile)
        stage = os.path.split(stageDir)[-1]
        stageNum = stage.replace('stage_', '').strip()
        outFile = '{}.tpr'.format(stage)
        tprFile = os.path.join(stageDir, outFile)
        if os.path.exists(tprFile): return tprFile
        groFile, topFile, _ = self.getPrevFinishedStageFiles(stage)


        if self.checkIfPrevTrj(stageNum):
            prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
        else:
            prevTrjStr = ''

        command = 'grompp -f %s -c %s -r %s -p ' \
                  '%s %s -o %s' % (os.path.abspath(mdpFile), groFile, groFile, topFile,
                                   prevTrjStr, outFile)
        #Manage warnings
        nWarns = self.countWarns(stageNum)
        print('{} warnings in stage {}'.format(nWarns, stageNum))
        if nWarns >= 1:
            command += ' -maxwarn {}'.format(nWarns)
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir)
        return tprFile

    def callMDRun(self, tprFile, saveTrj=True):
        stageDir = os.path.dirname(tprFile)
        stage = os.path.split(stageDir)[-1]
        gpuStr = ''
        if getattr(self, params.USE_GPU):
            gpuList = getattr(self, params.GPU_LIST).get().replace(' ', '')
            gpuStr = ' -gpu_id {}'.format(gpuList)

        command = 'mdrun -v -deffnm {}{} -nt {} -pin on -cpi -cpt {}'.format(stage, gpuStr, self.numberOfThreads.get(),
                                                                             self.cptTime.get())

        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir)
        trjFile = os.path.join(stageDir, '{}.trr'.format(stage))
        if os.path.exists(trjFile) and not saveTrj:
            os.remove(trjFile)

    def getPrevFinishedStageFiles(self, stage=None, reverse=False):
        '''Return the previous .gro and topology files if number stage is provided.
        If not, returns the ones of the lastest stage'''
        topFile = self.gromacsSystem.get().getTopologyFile()
        if stage:
            stageNum = stage.replace('stage_', '').strip()
            if stageNum == '1':
                groFile = os.path.abspath(self.gromacsSystem.get().getSystemFile())
                tprFile = None

            else:
                prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
                for file in os.listdir(prevDir):
                    if '.gro' in file:
                        groFile = os.path.join(prevDir, file)
                    elif '.tpr' in file:
                      tprFile = os.path.join(prevDir, file)
        else:
            stageDirs = natural_sort(glob.glob(self._getExtraPath('stage_*')), rev=reverse)
            for file in os.listdir(stageDirs[-1]):
                if '.gro' in file:
                    groFile = os.path.join(stageDirs[-1], file)
                elif '.tpr' in file:
                    tprFile = os.path.join(stageDirs[-1], file)

        return os.path.abspath(groFile),  os.path.abspath(topFile), tprFile

    def checkIfPrevTrj(self, stageNum):
        if stageNum == '1':
            return False
        else:
            prevDir = self._getExtraPath('stage_{}'.format(int(stageNum) - 1))
            for file in os.listdir(prevDir):
                if '.cpt' in file:
                  return os.path.join(prevDir, file)
        return False

    def getTrjFiles(self):
        trjFiles = []
        stagesDirs = natural_sort(glob.glob(self._getExtraPath('stage_*')), rev=True)
        for sDir in stagesDirs:
            cont = False
            for file in os.listdir(sDir):
                if '.trr' in file:
                  trjFiles.append(os.path.abspath(os.path.join(sDir, file)))
                  cont = True
            if not cont:
                break

        #Add previous trajectory if all stages in this protocol saved their trajectory (continuity)
        if len(trjFiles) == len(stagesDirs) and self.gromacsSystem.get().hasTrajectory() and self.prevTrj.get():
            prevTrjFile = os.path.abspath(self.gromacsSystem.get().getTrajectoryFile())
            conTrjFile = os.path.abspath(self._getTmpPath('previousTrajectory.trr'))

            command = 'trjconv -f {} -o {}'.format(prevTrjFile, conTrjFile)
            gromacsPlugin.runGromacs(self, 'gmx', command, cwd=self._getPath())

            trjFiles.append(conTrjFile)

        trjFiles.reverse()
        return trjFiles

    def concatTrjFiles(self, outTrj, tprFile):
        trjFiles = self.getTrjFiles()
        if len(trjFiles) > 0:
            tmpTrj = os.path.abspath(self._getTmpPath('concatenated.xtc'))
            #Concatenates trajectory
            command = 'trjcat -f {} -settime -o {} -cat'.format(' '.join(trjFiles), tmpTrj)
            gromacsPlugin.runGromacsPrintf(printfValues=['c'] * len(trjFiles),
                                           args=command, cwd=self._getPath())
            #Fixes and center trajectory
            command = 'trjconv -s {} -f {} -center -ur compact -pbc mol -o {}'.\
              format(os.path.abspath(tprFile), tmpTrj, outTrj)
            gromacsPlugin.runGromacsPrintf(printfValues=['Protein', 'System'] * len(trjFiles),
                                           args=command, cwd=self._getPath())
            return os.path.abspath(self._getPath(outTrj))
        return None

    def countWarns(self, stageNum):
        nWarns = 0
        for warn in self._warnings():
            if warn.split()[1] in ['all', str(stageNum)]:
                nWarns += 1
        return nWarns

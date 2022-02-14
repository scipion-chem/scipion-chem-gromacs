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
from pyworkflow.protocol.params import (LEVEL_ADVANCED, LEVEL_NORMAL, USE_GPU)
from pwem.protocols import EMProtocol

from pwchem.utils import natural_sort

from gromacs.objects import *
from gromacs.constants import *
from gromacs import Plugin as gromacsPlugin


class GromacsMDSimulation(EMProtocol):
    """
    This protocol will perform energy minimization on the system previosly prepared by the protocol "system prepartion".
    This step is necessary to energy minize the system in order to avoid unwanted conformations.
    """
    _label = 'Molecular dynamics simulation'
    _ensemTypes = ['Energy min', 'NVT',  'NPT']
    INTEGRATOR_OPTIONS = ['steep', 'cg', 'md']

    _thermostats = ['V-rescale']
    _barostats = ['Parrinello-Rahman']
    _coupleStyle = ['isotropic', 'semi-isotropic', 'anisotropic', 'constant area'] #check
    _restrainTypes = ['None', 'Ligand', 'Protein', 'Solute_heavy_atom', 'Solute']

    _paramNames = ['simTime', 'timeStep', 'saveTrj', 'trajInterval',
                   'temperature', 'tempRelaxCons', 'pressure', 'presRelaxCons', 'restrainForce']
    _enumParamNames = ['integrator', 'ensemType', 'thermostat', 'barostat', 'coupleStyle', 'restrains']
    _defParams = {'simTime': 100, 'timeStep': 0.002, 'saveTrj': False, 'trajInterval': 1.0,
                  'temperature': 300.0, 'tempRelaxCons': 0.1, 'integrator': 'md',
                  'pressure': 1.0, 'presRelaxCons': 2.0, 'restrainForce': 50.0,
                  'ensemType': 'NVT', 'thermostat': 'V-rescale', 'barostat': 'Parrinello-Rahman',
                  'coupleStyle': 'isotropic', 'restrains': 'None'}

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
      EMProtocol.__init__(self, **kwargs)


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):

        """ Define the input parameters that will be used.
        """

        form.addSection(label=Message.LABEL_INPUT)
        form.addHidden(params.USE_GPU, params.BooleanParam,
                      label='Use GPU for execution',
                      default=False,
                      help="GPU may have several cores. Set it one if "
                           "you don't know what we are talking about but you have a GPU."
                           "For DARC, first core index is 1, second 2, and so on. Write 0 if you do not want"
                           "to use GPU")

        form.addHidden(params.GPU_LIST, params.StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")

        form.addParam('gromacsSystem', params.PointerParam, label="Input Gromacs System: ",
                      pointerClass='GromacsSystem',
                      allowsNull=True,
                      help='Gromacs solvated system to be simulated')
        form.addParam('integrator', params.EnumParam,
                       label='Simulation integrator: ',
                       choices=self.INTEGRATOR_OPTIONS, default=0,
                       help='Type of integrator to use in simulation.')

        form.addParam('paramsFromFile', params.BooleanParam, default=False,
                      label="Use parameters from file (.mdp): ", expertLevel=LEVEL_ADVANCED,
                      help='Load simulation parameters from a local file.')
        form.addParam('ParamsFile', params.PathParam,
                      label='Parameters file: ',
                      condition="paramsFromFile",
                      allowsnull=False,
                      help="Input parameters file which will be used to run energy minimization")

        group = form.addGroup('Trajectory')
        group.addParam('saveTrj', params.BooleanParam, default=self._defParams['saveTrj'],
                       label="Save trajectory: ",
                       help='Save trajectory of the atoms during stage simulation.'
                            'The output will concatenate those trajectories which appear after the last stage '
                            'where the trajectory was not saved.')
        group.addParam('trajInterval', params.FloatParam, default=self._defParams['trajInterval'],
                       label='Interval time (ps):', condition='saveTrj',
                       help='Time between each frame recorded in the simulation (ps)')

        group = form.addGroup('Simulation time')
        group.addParam('simTime', params.FloatParam, default=100,
                       label='Simulation time (ps):',
                       help='Total time of the simulation stage (ps)')
        group.addParam('timeStep', params.FloatParam, default=0.002,
                       label='Simulation time steps, dt (ps):', expertLevel=LEVEL_ADVANCED,
                       help='Time of the steps for simulation, dt (ps)')

        group = form.addGroup('Ensemble')
        group.addParam('ensemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
                       help='Type of simulation to perform in the step: Energy minimization, NVT or NPT')

        line = group.addLine('Temperature settings: ', condition='ensemType!=0',
                             help='Temperature during the simulation (K)\nThermostat type\n'
                                  'Relaxation time constant for thermostat (ps)')
        line.addParam('temperature', params.FloatParam, default=300, condition='ensemType!=0',
                      label='Temperature: ')
        line.addParam('thermostat', params.EnumParam, default=0, condition='ensemType!=0',
                      label='Thermostat type: ', choices=self._thermostats, expertLevel=LEVEL_ADVANCED)
        line.addParam('tempRelaxCons', params.FloatParam, default=0.1,
                      label='Temperature time constant (ps): ', expertLevel=LEVEL_ADVANCED)

        line = group.addLine('Pressure settings: ', condition='ensemType==2',
                             help='Pressure during the simulation (bar)\nBarostat type\n'
                                  'Relaxation time constant for barostat (ps)')
        line.addParam('pressure', params.FloatParam, default=1.0,
                      label='   Pressure (bar):   ')
        line.addParam('barostat', params.EnumParam, default=0,
                      label='  Barostat type:   ', choices=self._barostats, expertLevel=LEVEL_ADVANCED)
        line.addParam('presRelaxCons', params.FloatParam, default=2.0,
                      label='   Pressure time constant (ps):   ', expertLevel=LEVEL_ADVANCED)
        group.addParam('coupleStyle', params.EnumParam, default=0, condition='ensemType==2',
                       label='Pressure coupling style: ', choices=self._coupleStyle,
                       expertLevel=LEVEL_ADVANCED)

        group = form.addGroup('Restrains')
        group.addParam('restrains', params.EnumParam, default=0,
                       label='Restrains: ', choices=self._restrainTypes,
                       help='Restrain movement of specific groups of atoms')
        group.addParam('restrainForce', params.FloatParam, default=50,
                       label='Restrain force constant: ', condition='restrains!=0',
                       help='Restrain force applied to the selection (kcal/mol/Å2)')

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
        i=1
        for wStep in self.workFlowSteps.get().split('\n'):
            if wStep:
              self._insertFunctionStep('simulateStageStep', wStep, i)
              i+=1
        self._insertFunctionStep('createOutputStep')

    def simulateStageStep(self, wStep, i):
      if wStep in ['', None]:
          msjDic = self.createMSJDic()
      else:
          msjDic = eval(wStep)
      mdpFile = self.generateMDPFile(msjDic, i)

      tprFile = self.callGROMPP(mdpFile)
      self.callMDRun(tprFile, saveTrj=msjDic['saveTrj'])

    def createOutputStep(self):
        lastGroFile, lastTopoFile, lastTprFile = self.getPrevFinishedStageFiles()
        localGroFile, localTopFile = self._getPath('outputSystem.gro'), self._getPath('systemTopology.top')
        shutil.copyfile(lastGroFile, localGroFile)
        shutil.copyfile(lastTopoFile, localTopFile)

        outTrj = self.concatTrjFiles(outTrj='outputTrajectory.xtc', tprFile=lastTprFile)

        outSystem = GromacsSystem(filename=localGroFile)
        outSystem.setTopologyFile(localTopFile)
        outSystem.setTrajectoryFile(outTrj)

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

    def createSummary(self, msjDic=None):
        '''Creates the displayed summary from the internal state of the steps'''
        sumStr = ''
        if not msjDic:
            for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
                if dicLine != '':
                    msjDic = eval(dicLine)
                    msjDic = self.addDefaultForMissing(msjDic)
                    method, ensemType = msjDic['thermostat'], msjDic['ensemType']
                    sumStr += '{}) Sim. time ({}): {} ps, {} ensemble'.\
                      format(i+1, msjDic['integrator'], msjDic['simTime'], ensemType)

                    if msjDic['restrains'] != 'None':
                        sumStr += ', restrain on {}'.format(msjDic['restrains'])
                    sumStr += ', {} K\n'.format(msjDic['temperature'])
        else:
            msjDic = self.addDefaultForMissing(msjDic)
            method, ensemType = msjDic['thermostat'], msjDic['ensemType']
            sumStr += 'Sim. time ({}): {} ps, {} ensemble'. \
              format(msjDic['simTime'], msjDic['integrator'], ensemType, method)

            if msjDic['restrains'] != 'None':
              sumStr += ', restrain on {}'.format(msjDic['restrains'])
            sumStr += ', {} K\n'.format(msjDic['temperature'])
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
        os.mkdir(stageDir)
        mdpFile = os.path.join(stageDir, 'stage_{}.mdp'.format(mdpStage))

        restr = msjDic['restrains']
        if restr != 'None':
            #TODO: generate new restriction itp file if non considered before, add to topology, add line
            restrStr = RESTR_STR
        else:
            restrStr = ''

        integStr = msjDic['integrator']
        tSteps = msjDic['timeStep']
        nSteps = round(msjDic['simTime'] / tSteps)
        if msjDic['ensemType'] == 'Energy min':
            tStepsStr = TSTEP_EM.format(tSteps)
            dispCorrStr, outControlStr = '', ''
            bondParStr, electroStr = '', ''
            tempStr, presStr = '', ''
            velStr, velParStr = 'no', ''
        else:
            tStepsStr = TSTEP_EQ.format(tSteps)
            dispCorrStr = DISP_CORR
            electroStr = ELECTROSTATICS
            tempStr = TEMP_SETTING.format(msjDic['thermostat'],
                                          *[msjDic['temperature']]*2, *[msjDic['tempRelaxCons']]*2)
            if msjDic['ensemType'] == 'NVT':
                presStr = ''
            elif msjDic['ensemType'] == 'NPT':
                presStr = PRES_SETTING.format(msjDic['barostat'], msjDic['coupleStyle'],
                                              msjDic['pressure'], msjDic['presRelaxCons'])

            if self.checkIfPrevTrj(mdpStage):
                bondParStr = BONDED_PARAMS.format('yes')
                velStr, velParStr = 'no', ''
            else:
                bondParStr = BONDED_PARAMS.format('no')
                velStr = 'yes'
                velParStr = VEL_GEN.format(msjDic['temperature'])

        nTrj = round(msjDic['trajInterval'] / tSteps)
        if msjDic['saveTrj']:
            outControlStr = OUTPUT_CONTROL.format(*[nTrj]*4)
        else:
            outControlStr = OUTPUT_CONTROL.format(*[0]*3, nTrj)

        title = 'Stage {}: {}, {} ps'.format(mdpStage, msjDic['ensemType'], msjDic['simTime'])
        mdpStr = MDP_STR.format(title, restrStr, integStr, nSteps, tStepsStr, dispCorrStr, outControlStr, bondParStr,
                                electroStr, tempStr, presStr, velStr, velParStr)

        with open(mdpFile, 'w') as f:
            f.write(mdpStr)

        return mdpFile

    def callGROMPP(self, mdpFile):
        stageDir = os.path.dirname(mdpFile)
        stage = os.path.split(stageDir)[-1]
        stageNum = stage.replace('stage_', '').strip()
        groFile, topFile, _ = self.getPrevFinishedStageFiles(stage)
        outFile = '{}.tpr'.format(stage)

        if self.checkIfPrevTrj(stageNum):
            prevTrjStr = '-t ' + os.path.abspath(self.checkIfPrevTrj(stageNum))
        else:
            prevTrjStr = ''

        command = 'grompp -f %s -c %s -r %s -p ' \
                  '%s %s -o %s' % (os.path.abspath(mdpFile), groFile, groFile, topFile,
                                   prevTrjStr, outFile)
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir)
        return os.path.join(stageDir, outFile)

    def callMDRun(self, tprFile, saveTrj=True):
        stageDir = os.path.dirname(tprFile)
        stage = os.path.split(stageDir)[-1]
        gpuStr = ''
        if getattr(self, params.USE_GPU):
            gpuStr = ' -nb gpu'

        command = 'mdrun -v -deffnm {} {}'.format(stage, gpuStr)
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir)
        if not saveTrj:
            trjFile = os.path.join(stageDir, '{}.trr'.format(stage))
            os.remove(trjFile)

    def getPrevFinishedStageFiles(self, stage=None, reverse=False):
        '''Return the previous .gro and topology files if number stage is provided.
        If not, returns the ones of the lastest stage'''
        topFile = self.gromacsSystem.get().getTopologyFile()
        #todo: copy topology to extra and modify there if new restrictions
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

        trjFiles.reverse()
        return trjFiles

    def concatTrjFiles(self, outTrj, tprFile):
        trjFiles = self.getTrjFiles()
        tmpTrj = os.path.abspath(self._getTmpPath('concatenated.xtc'))
        #Concatenates trajectory
        command = 'trjcat -f {} -settime -o {}'.format(' '.join(trjFiles), tmpTrj)
        gromacsPlugin.runGromacsPrintf(self, 'gmx', printfValues=['c'] * len(trjFiles),
                                       args=command, cwd=self._getPath())
        #Fixes and center trajectory
        command = 'trjconv -s {} -f {} -center -ur compact -pbc mol -o {}'.\
          format(os.path.abspath(tprFile), tmpTrj, outTrj)
        gromacsPlugin.runGromacsPrintf(self, 'gmx', printfValues=['Protein\nSystem\n'] * len(trjFiles),
                                       args=command, cwd=self._getPath())
        return self._getPath(outTrj)

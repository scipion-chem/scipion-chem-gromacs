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
import glob, uuid
import os.path

from pyworkflow.object import String
from pyworkflow.protocol import params
from pyworkflow.utils import Message, runJob, createLink
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pwchem.utils import natural_sort
from pwem.convert import AtomicStructHandler

from gromacs.objects import *
from gromacs.constants import *
from gromacs import Plugin as gromacsPlugin, GromacsSystem

from multiprocessing import cpu_count

class GromacsMDSimulation(EMProtocol):
    """
    This protocol will perform energy minimization on the system previosly prepared by the protocol "system prepartion".

    AI Generated:

        GromacsMDSimulation

        Overview
        --------
        This protocol performs molecular dynamics (MD) simulations using GROMACS
        within the Scipion-Chem framework.

        It supports energy minimization, NVT, and NPT simulations, allowing flexible
        configuration of thermodynamic conditions, restraints, and simulation stages.

        The protocol enables multi-step workflows, where each stage can be customized
        and executed sequentially, producing trajectories and system states for
        downstream analysis.

        Inputs
        ------
        gromacsSystem:
            A previously prepared GROMACS system containing:
            - Coordinates (GRO/PDB)
            - Topology (TOP)
            - Optional trajectory and index files

        prevTrj:
            Whether to concatenate the previous trajectory from the input system

        cptTime:
            Time interval (minutes) for writing checkpoint files

        gmxMPI:
            Whether to use MPI execution for GROMACS

        Workflow
        --------
        1. Workflow definition
           - Reads user-defined simulation stages from workflow steps
           - Each stage is represented as a parameter dictionary
           - Supports sequential multi-stage simulations

        2. Stage preparation
           - Generates MDP configuration files for each stage
           - Configures:
             - Integrator (steep, cg, md)
             - Temperature and thermostat
             - Pressure and barostat (for NPT)
             - Time step and simulation length
             - Output frequency and trajectory saving

        3. Restraints handling
           - Supports:
             - Structural ROI restraints
             - Residue-based restraints
             - Custom index groups (make_ndx)
           - Generates index files and applies force constants

        4. Preprocessing (GROMPP)
           - Converts MDP + structure + topology into TPR files
           - Handles warnings and allows controlled tolerance

        5. Simulation execution (MDRUN)
           - Runs MD simulation using CPU or GPU
           - Supports MPI parallelization
           - Generates:
             - Trajectory files (.trr/.xtc)
             - Checkpoints (.cpt)
             - Updated structures (.gro)

        6. Trajectory management
           - Concatenates trajectories across stages
           - Optionally merges with previous system trajectory
           - Applies centering and periodic boundary corrections

        7. Output generation
           - Copies final structure and topology
           - Builds final GromacsSystem object
           - Attaches trajectory and index files if available

        Output
        ------
        outputSystem:
            GromacsSystem containing:
            - Final coordinates (GRO)
            - Topology (TOP)
            - Simulation parameters (TPR)
            - Optional trajectory (XTC)
            - Index and restraint definitions

        Summary
        -------
        This protocol provides a complete framework for molecular dynamics simulations,
        enabling:
        - multi-stage simulation workflows (minimization, equilibration, production)
        - flexible thermodynamic control (NVT/NPT)
        - GPU and MPI acceleration
        - trajectory generation and concatenation
        - application of structural and custom restraints
        - reproducible and modular simulation pipelines

        Notes
        -----
        - Supports energy minimization and MD simulations in a unified workflow
        - Compatible with GPU acceleration and MPI execution
        - Requires a properly prepared GROMACS system as input
        - Large trajectories may require significant storage
        - Includes validation and warnings for simulation stability and parameter consistency

    """
    _label = 'Run MD simulation'
    _ensemTypes = ['Energy min', 'NVT',  'NPT']

    _integrators = ['steep', 'cg']
    _thermostats = ['no', 'Berendsen', 'Nose-Hoover', 'Andersen', 'Andersen-massive', 'V-rescale']
    _barostats = ['no', 'Berendsen', 'Parrinello-Rahman', 'C-rescale']
    _coupleStyle = ['isotropic', 'semiisotropic']
    _restraints = ['Structural ROI', 'Residues', 'Custom make_ndx command']

    _omitParamNames = ['useGpu', 'gpuList', 'gromacsSystem', 'restrainROIs',
                       'runName', 'runMode', 'insertStep', 'summarySteps', 'deleteStep', 'watchStep',
                       'workFlowSteps', 'hostName', 'numberOfThreads', 'numberOfMpi']

    # -------------------------- DEFINE constants ----------------------------
    def __init__(self, **kwargs):
      EMProtocol.__init__(self, **kwargs)

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        cpus = cpu_count()//2 # don't use everything
        form.addParallelSection(threads=cpus, mpi=1)

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

        form.addParam('gmxMPI', params.BooleanParam, default=False,
                       label="Use MPI program: ",
                       help='Use MPI program during simulation stage.')

        group = form.addGroup('Ensemble')
        group.addParam('ensemType', params.EnumParam,
                       label='Simulation type: ',
                       choices=self._ensemTypes, default=0,
                       help='Type of simulation to perform in the step: Energy minimization, NVT or NPT\n'
                            'https://manual.gromacs.org/5.1.1/user-guide/mdp-options.html')

        group.addParam('integrator', params.EnumParam, label='Simulation integrator: ', condition='ensemType==0',
                      choices=self._integrators, default=0, help='Type of integrator to use in simulation.')

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
        line.addParam('barostat', params.EnumParam, default=3,
                      label='  Barostat type:   ', choices=self._barostats)
        line.addParam('presRelaxCons', params.FloatParam, default=2.0,
                      label='   Pressure constant (ps)[tau-p]:   ', expertLevel=params.LEVEL_ADVANCED)
        line.addParam('presCouple', params.IntParam, default=-1,
                      label='Coupling frequency [nstpcouple]: ', expertLevel=params.LEVEL_ADVANCED)
        group.addParam('coupleStyle', params.EnumParam, default=0, condition='ensemType==2',
                      label='Pressure coupling style: ', choices=self._coupleStyle,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Semiisotropic is recomemded for membrane proteins')

        group = form.addGroup('Trajectory', condition='ensemType!=0')
        group.addParam('saveTrj', params.BooleanParam, default=False,
                       label="Save trajectory: ", condition='ensemType!=0',
                       help='Save trajectory of the atoms during simulation stage.'
                            'The output will concatenate those trajectories which appear after the last stage '
                            'where the trajectory was not saved.')
        group.addParam('trajInterval', params.FloatParam, default=1.0,
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
                       label='Initial step-size (nm)[emstep]: ', condition='ensemType==0',
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
        group.addParam('restraintOptions', params.EnumParam, label='Create new restraints group from: ',
                       choices=self._restraints, default=0, expertLevel=params.LEVEL_ADVANCED,
                       help='Type of restraint group to create')
        
        group.addParam('restrainROIs', params.PointerParam, label='Input set of ROIs: ', pointerClass='SetOfStructROIs', 
                       condition='restraintOptions==0', expertLevel=params.LEVEL_ADVANCED, allowsNull=True,
                       help='Restraint movement of the groups of atoms included in the specified structural ROI.')
        group.addParam('restrainROI', params.StringParam, label='ROI to retrain: ',
                       condition='restraintOptions==0', expertLevel=params.LEVEL_ADVANCED,
                       help='Restraint movement of the groups of atoms included in the specified structural ROI.')
        group.addParam('restraintROIInfo', params.LabelParam, label='Add ROI restraint: ',
                       condition='restraintOptions == 0', expertLevel=params.LEVEL_ADVANCED,
                       help='Generates the specified index group from the atoms in contact with the selected ROI')

        group.addParam('restrainChain', params.StringParam, label='Restrain in structure chain: ',
                       condition='restraintOptions == 1', expertLevel=params.LEVEL_ADVANCED,
                       help='Select the chain where the residues to restrain are')
        group.addParam('restrainResidue', params.StringParam, label='Restrain chain residues : ',
                       condition='restraintOptions == 1', expertLevel=params.LEVEL_ADVANCED,
                       help='Restraint movement of the groups of atoms included in the selected residues.')
        group.addParam('restraintResidueInfo', params.LabelParam, label='Add residue restraint: ',
                       condition='restraintOptions == 1', expertLevel=params.LEVEL_ADVANCED,
                       help='Generates the specified index group from the atoms in the selected residues')

        group.addParam('restraintCommand', params.StringParam, default='', label='Enter custom index command: ',
                       expertLevel=params.LEVEL_ADVANCED, condition='restraintOptions==2',
                       help='Define custom atom groups using raw GROMACS make_ndx syntax.\n\n'
                            'To execute, type your command here and click this wizard button. '
                            'The new group will be generated and can be verified/selected in the '
                            '"Choose restraints group" wizard below.\n\n'
                            'Syntax Cheat Sheet:\n'
                            '  & (AND), | (OR), ! (NOT)\n'
                            '  nr / "name"  -> Select by group number or exact name in quotes (e.g., 1 or "chainA")\n'
                            '  a [name]     -> Select by atom name (e.g., a C*)\n'
                            '  r [name/nr]  -> Select by residue (e.g., r 1-15 or r LYS)\n'
                            '  chain [char] -> Select by chain identifier (e.g., chain A)\n\n'
                            'Examples:\n'
                            '  1 & ! a H* (Protein group without Hydrogens)\n'
                            '  "chainA" & ! a H* (Chain A without Hydrogens)\n'
                            '  17 & t C            (Only Carbon atom types inside chainA)\n'
                            '  1 | 13              (Combines Protein and SOL groups)')
        
        group.addParam('restraints', params.StringParam, default='None', label='Choose restraints group: ',
                       help='Restraint movement of specific groups of atoms. You can check the existing groups of '
                            'your system using the wizard or even create new groups in the advanced level using '
                            'make_ndx gromacs commands.')
        group.addParam('restraintForce', params.FloatParam, default=50, label='Restraint force constant: ',
                       help='Restraint force applied to the selection (kcal/mol/Å2)')

        group = form.addGroup('Summary')
        group.addParam('insertStep', params.StringParam, default='',
                       label='Insert step number: ',
                       help='Insert the defined step into the workflow on the defined position.\n'
                            'The default (when empty) is the last position')
        group.addParam('summarySteps', params.TextParam, width=100, readOnly=True,
                       label='Summary of steps',
                       help='Summary of the defined steps. \nManual modification will have no '
                            'effect, use the wizards to add / delete the steps')
        group.addParam('deleteStep', params.StringParam, default='',
                       label='Delete step number: ',
                       help='Delete the step of the specified index from the workflow.')
        group.addParam('watchStep', params.StringParam, default='',
                       label='Watch parameters of step number: ',
                       help='Watch the parameters step of the specified index from the workflow...\n'
                            'This might be useful if you want to change some parameters of a predefined step. '
                            'However, the parameters are not changed until you add the new step (and probably '
                            'you may want to delete the previous unchanged step)')
        group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')

        form.addSection(label='Input Pointers')
        form.addParam('inputPointerLabels', params.LabelParam, important=True,
                      label='Records of inputs. Do not modify manually',
                      help='This is a list of the input pointer to keep track of the inputs received.\n'
                           'It is automatically updated with the first section wizards.\n'
                           'Manual modification (adding inputs from the lens) will have no actual impact on the '
                           'protocol performance')
        form.addParam('inputPointers', params.MultiPointerParam, pointerClass="Sequence, AtomStruct",
                      label='Input Pointers: ', allowsNull=True,
                      help='This is a list of the input pointer to keep track of the inputs received.\n'
                           'It is automatically updated with the first section wizards.\n'
                           'Manual modification (adding inputs from the lens) will have no actual impact on the '
                           'protocol performance')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self.createGUISummary()
        i = 1
        for wStep in self.workFlowSteps.get().strip().split('\n'):
            self._insertFunctionStep(self.simulateStageStep, wStep, i)
            i += 1
        self._insertFunctionStep(self.createOutputStep)

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
        shutil.copyfile(lastGroFile, localGroFile), shutil.copyfile(lastTopoFile, localTopFile)
        outTrj = self.concatTrjFiles(outTrj='outputTrajectory.xtc', tprFile=lastTprFile)
        localPdbFile = self._getPath('outputSystem.pdb')
        self._convertGroToPdbNoWat(localGroFile, localPdbFile)
        finalAtomStruct = AtomStruct(filename=os.path.relpath(localPdbFile))

        outSystem = GromacsSystem(filename=localGroFile, oriStructFile=oriGroFile, tprFile=lastTprFile)
        outSystem.setTopologyFile(localTopFile)
        outSystem.setLigTopologyFile(self.gromacsSystem.get().getLigTopologyFile())
        outSystem.setLigandID(self.gromacsSystem.get().getLigandID())
        outSystem.setChainNames(self.gromacsSystem.get().getChainNames())
        outSystem.setChainLengths(self.gromacsSystem.get().getChainLengths())
        if outTrj:
            outSystem.setTrajectoryFile(outTrj)
            outSystem.readTrjInfo(protocol=self, outDir=self._getExtraPath())

        # set and clean indexFiles
        customIndex = gromacsPlugin.getCustomIndexFile(self)
        indexFile = self._getExtraPath('indexes.ndx')
        if os.path.exists(customIndex):
            shutil.copy(customIndex, indexFile)
        else:
            shutil.copy(self.gromacsSystem.get().getIndexFile(), indexFile)
        self.cleanCustomIndex()
        outSystem.setIndexFile(indexFile)

        self._defineOutputs(outputSystem=outSystem, lastFrameStruct=finalAtomStruct)


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
  
        if msjDic['restraints'] not in ['', 'None']:
          sumStr += ', restraint on {}'.format(msjDic['restraints'])
        sumStr += ', {} K\n'.format(msjDic['temperature'])
        return sumStr

    def createSummary(self, msjDic=None):
        '''Creates the displayed summary from the internal state of the steps'''
        if not msjDic:
            sumStr = ''
            for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
                dicLine = dicLine.strip()
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

            if msjDic['ensemType'] != 'Energy min':
                if msjDic['thermostat'] not in ['Berendsen', 'V-rescale'] and not prevTrj and msjDic['saveTrj']:
                    warns.append(f'\nStep {step+1} : Berendsen and V-rescale are the thermostat recommended for '
                                 f'system equilibration, {msjDic["thermostat"]} might not be the best option for '
                                 f'the first trajectory saved')

                if prevTrj and msjDic['thermostat'] == 'Berendsen':
                    warns.append(f'\nStep {step+1} : Berendsen is the thermostat recommended for system equilibration '
                                 f'only, it might not be the best option for later trajectories saved')
                tCoup = msjDic['timeNeigh'] if msjDic['tempCouple'] == -1 else msjDic['tempCouple']
                if msjDic['thermostat'] == 'Nose-Hoover' and 20*tCoup*msjDic['timeStep'] > msjDic['tempRelaxCons']:
                    warns.append(f'\nStep {step+1} : For proper integration of the Nose-Hoover thermostat, tau-t '
                                 f'({msjDic["tempRelaxCons"]}) should '
                                 f'be at least 20 times larger than nsttcouple*dt ({tCoup}*{msjDic["timeStep"]})')

                if msjDic['ensemType'] == 'NPT':
                  if msjDic['barostat'] not in ['Berendsen', 'C-rescale'] and not prevTrj and msjDic['saveTrj']:
                    warns.append(f'\nStep {step + 1} : Berendsen and C-rescale are the barostats recommended for '
                                 f'system equilibration, {msjDic["barostat"]} might not be the best option for '
                                 f'the first trajectory saved')
                  if prevTrj and msjDic['barostat'] == 'Berendsen':
                    warns.append(f'\nStep {step + 1} : Berendsen is the barostat recommended for system equilibration '
                                 f'only, it might not be the best option for later trajectories saved')

            if msjDic['saveTrj']:
                prevTrj = True

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
                  vals.append(f'Step {step+1} : Andersen temperature control not supported for integrator md.')
        return vals

######################## UTILS ##################################

    def getPrevPointersIds(self, prevPointers):
      ids = []
      for p in prevPointers:
        ids.append(p.get().getObjId())
      return ids

    def _convertGroToPdb(self, groFile, pdbFile):
        """ Helper function to convert GRO to PDB using GROMACS editconf """
        args = f'-f {groFile} -o {pdbFile}'
        gromacsPlugin.runGromacs(self, 'gmx editconf', args)

    def _convertGroToPdbNoWat(self, groFile, pdbFile):
        """ Helper function to convert GRO to PDB while preserving chain labels """
        inpSystem = self.gromacsSystem.get()
        modelChains = inpSystem.getChainNames()
        indexFile = inpSystem.getIndexFile()

        # Case 1: Single chain processing
        if len(modelChains) == 1:
            printGroup = ['Protein']
            if inpSystem.hasLig():
                printGroup = [f'Protein_{inpSystem.getLigandID()}']

            self._runEditconf(groFile, indexFile, pdbFile, modelChains[0], printGroup)
            return

        # Case 2: Multiple chains - process each separately and combine
        allPdbs = []

        # Extract each protein chain with proper chain label
        for chainId in modelChains:
            chainPdb = self._getTmpPath(f'chain_{chainId}.pdb')
            allPdbs.append(chainPdb)
            self._runEditconf(groFile, indexFile, chainPdb, chainId, [f'chain{chainId}'])

        # Add ligand if present
        if inpSystem.hasLig():
            ligPdb = self._getTmpPath('ligand.pdb')
            allPdbs.append(ligPdb)
            self._runEditconf(groFile, indexFile, ligPdb, 'L', [inpSystem.getLigandID()])

        # Combine all individual PDB pieces into the final file
        self._combinePdbFiles(allPdbs, pdbFile)

    def _runEditconf(self, groFile, indexFile, outFile, label, printGroup):
        """ Runs the Gromacs editconf command with specified label and group """
        params = " editconf -f {} -n {} -o {} -label {}".format(
            os.path.abspath(groFile),
            os.path.abspath(indexFile),
            os.path.abspath(outFile),
            label)
        gromacsPlugin.runGromacsPrintf(self, printfValues=printGroup,
                                       args=params, cwd=self._getPath())

    def _combinePdbFiles(self, pdbFiles, targetPdbFile):
        """ Merges structural data lines from multiple PDB files into one """
        with open(targetPdbFile, 'w') as outFile:
            firstFile = True

            for pdb in pdbFiles:
                with open(pdb, 'r') as inFile:
                    for line in inFile:
                        # Keep the unit cell info (CRYST1) only from the very first file
                        if firstFile and line.startswith('CRYST1'):
                            outFile.write(line)
                            firstFile = False

                        # Only keep core coordinate records
                        if line.startswith(('ATOM', 'HETATM')):
                            outFile.write(line)

            outFile.write('END\n')

    def countSteps(self):
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

    def getStageParamsDic(self, type='All'):
      '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
      Type'''
      paramsDic = {}
      for paramName, param in self._definition.iterAllParams():
        if paramName not in self._omitParamNames and not isinstance(param, params.Group) and not isinstance(param, params.Line):
          if type == 'All':
            paramsDic[paramName] = param
          elif type == 'Enum' and isinstance(param, params.EnumParam):
            paramsDic[paramName] = param
          elif type == 'Normal' and not isinstance(param, params.EnumParam):
            paramsDic[paramName] = param
      return paramsDic

    def createMSJDic(self):
        msjDic = {}
        for pName in self.getStageParamsDic(type='Normal').keys():
            if hasattr(self, pName):
                msjDic[pName] = getattr(self, pName).get()
            else:
                print('Something is wrong with parameter ', pName)

        for pName in self.getStageParamsDic(type='Enum').keys():
            if hasattr(self, pName):
                msjDic[pName] = self.getEnumText(pName)
            else:
                print('Something is wrong with parameter ', pName)
        return msjDic

    def addDefaultForMissing(self, msjDic):
      '''Add default values for missing parameters in the msjDic'''
      paramDic = self.getStageParamsDic()
      for pName in paramDic.keys():
        if pName not in msjDic:
          msjDic[pName] = paramDic[pName].default
      return msjDic

    def generateMDPFile(self, msjDic, mdpStage):
        stageDir = self._getExtraPath('stage_{}'.format(mdpStage))
        if not os.path.exists(stageDir):
            os.mkdir(stageDir)

        mdpFile = os.path.join(stageDir, 'stage_{}.mdp'.format(mdpStage))
        if os.path.exists(mdpFile): return mdpFile

        indexFile = os.path.abspath(gromacsPlugin.ensureIndexFile(self))

        if msjDic['restraints'].strip() not in ('', 'None'):
            rSuffix = f"{msjDic['restraints']}_stg{mdpStage}"
            groupNr = gromacsPlugin.translateNamesToIndexGroup(self, [msjDic['restraints']])

            newSuffixes = self.gromacsSystem.get().\
              defineNewRestriction(self, index=groupNr, energy=msjDic['restraintForce'], restraintSuffix=rSuffix,
                                   outDir=stageDir, indexFile=indexFile)

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
                if msjDic['coupleStyle'] == 'isotropic':
                    presStr = PRES_SETTING.format(msjDic['barostat'], msjDic['coupleStyle'],
                                              msjDic['pressure'], msjDic['presRelaxCons'],
                                              msjDic['presCouple'])
                else:
                    presStr = PRES_SETTING_SEMI.format(msjDic['barostat'], msjDic['coupleStyle'],
                                              msjDic['pressure'], msjDic['pressure'], msjDic['presRelaxCons'],
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

        localTop = os.path.join(stageDir, os.path.split(topFile)[-1])
        if not os.path.exists(localTop):
          os.link(topFile, localTop)

        command = 'grompp -f %s -c %s -r %s -p ' \
                  '%s %s -o %s' % (os.path.abspath(mdpFile), groFile, groFile, os.path.split(topFile)[-1],
                                   prevTrjStr, outFile)

        ligTopFile = self.gromacsSystem.get().getLigTopologyFile()
        if ligTopFile:
          lTopFile = os.path.join(stageDir, os.path.split(ligTopFile)[-1])
          os.link(ligTopFile, lTopFile)

        #Manage warnings
        nWarns = self.countWarns(stageNum)
        print('{} warnings in stage {}'.format(nWarns, stageNum))
        if nWarns >= 1:
            command += ' -maxwarn {}'.format(nWarns)
        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir, numberOfMpi=0)
        return tprFile

    def callMDRun(self, tprFile, saveTrj=True):
        stageDir = os.path.dirname(tprFile)
        stage = os.path.split(stageDir)[-1]
        if getattr(self, params.USE_GPU):
            gpuList = getattr(self, params.GPU_LIST).get().replace(' ', '')
            gpuStr = f' -nb gpu -gpu_id {gpuList}'
        else:
            gpuStr = ' -nb cpu'

        gmxMPIStr = f'-ntomp {self.numberOfThreads.get()}' if self.gmxMPI.get() else f'-nt {self.numberOfThreads.get()}'
        command = f'mdrun -v -deffnm {stage}{gpuStr} {gmxMPIStr} -pin on -cpi -cpt {self.cptTime.get()}'

        gromacsPlugin.runGromacs(self, 'gmx', command, cwd=stageDir, mpi=self.gmxMPI.get())
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
            gromacsPlugin.runGromacsPrintf(self, printfValues=['c'] * len(trjFiles),
                                           args=command, cwd=self._getPath())
            #Fixes and center trajectory
            command = 'trjconv -s {} -f {} -center -ur compact -pbc mol -o {}'.\
              format(os.path.abspath(tprFile), tmpTrj, outTrj)
            gromacsPlugin.runGromacsPrintf(self, printfValues=['Protein', 'System'] * len(trjFiles),
                                           args=command, cwd=self._getPath())
            return os.path.abspath(self._getPath(outTrj))
        return None

    def countWarns(self, stageNum):
        nWarns = 0
        for warn in self._warnings():
            if warn.split()[1] in ['all', str(stageNum)]:
                nWarns += 1
        return nWarns

    def cleanCustomIndex(self):
        tmpPath = self.getProject().getTmpPath()
        for customInxFile in glob.iglob(os.path.join(tmpPath, "*custom_indexes.ndx*")):
            if os.path.isfile(customInxFile):
                os.remove(customInxFile)
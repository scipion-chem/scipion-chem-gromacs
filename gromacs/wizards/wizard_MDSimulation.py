# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

"""
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
import json, os
from pyworkflow.gui import ListTreeProviderString, dialog
import pyworkflow.object as pwobj

from pwem.objects import Pointer, String

from pwchem.utils import groupConsecutiveIdxs
from pwchem.wizards import AddElementSummaryWizard, DeleteElementWizard, VariableWizard, SelectElementWizard, \
    WatchElementWizard

from ..protocols import GromacsSystemPrep, GromacsMDSimulation
from gromacs import Plugin as gromacsPlugin

SelectElementWizard().addTarget(protocol=GromacsSystemPrep,
                                targets=['inputLigand'],
                                inputs=['inputSetOfMols'],
                                outputs=['inputLigand'])

AddElementSummaryWizard().addTarget(protocol=GromacsMDSimulation,
                             targets=['insertStep'],
                             inputs=['insertStep'],
                             outputs=['workFlowSteps', 'summarySteps'])

DeleteElementWizard().addTarget(protocol=GromacsMDSimulation,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])

WatchElementWizard().addTarget(protocol=GromacsMDSimulation,
                                targets=['watchStep'],
                                inputs=['watchStep'],
                                outputs=['workFlowSteps', 'summarySteps'])

class GromacsCheckIndexWizard(VariableWizard):
    """Watch the groups contained in the input Gromacs system"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        _, outputParam = self.getInputOutput(form)

        indexFile = gromacsPlugin.ensureIndexFile(protocol)
        groups = gromacsPlugin.parseIndexFile(protocol, indexFile)

        finalList = [String('-1: None')]
        for index, name in groups.items():
            finalList.append(String('{}: {}'.format(index, name)))
        provider = ListTreeProviderString(finalList)
        dlg = dialog.ListDialog(form.root, "System groups", provider,
                                "Select one of the groups")
        form.setVar(outputParam[0], dlg.values[0].get().split(': ')[1])

GromacsCheckIndexWizard().addTarget(protocol=GromacsMDSimulation,
                               targets=['restraints'],
                               inputs=['gromacsSystem'],
                               outputs=['restraints'])

class GromacsCustomIndexWizard(GromacsCheckIndexWizard):
    """Watch the groups contained in the input Gromacs system and allows the modification"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParam, _ = self.getInputOutput(form)
        inpSystem = getattr(protocol, inputParam[0]).get()

        inIndex, outIndex = gromacsPlugin.ensureIndexFile(protocol), gromacsPlugin.getCustomIndexFile(protocol)

        inCommand = getattr(protocol, inputParam[1]).get()
        inCommand = ' | '.join(map(str, gromacsPlugin.translateNamesToIndexGroup(protocol, inCommand.split())))

        gromacsPlugin.createIndexFile(protocol, inpSystem, inputCommands=[inCommand], inIndex=inIndex, outIndex=outIndex)

GromacsCustomIndexWizard().addTarget(protocol=GromacsMDSimulation,
                               targets=['restraintCommand'],
                               inputs=['gromacsSystem', 'restraintCommand'],
                               outputs=['restraints'])

SelectElementWizard().addTarget(protocol=GromacsMDSimulation,
                                targets=['restrainROI'],
                                inputs=['restrainROIs'],
                                outputs=['restrainROI'])


class SelectChainWizardGromacs(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        try:
          inputObj = getattr(protocol, inputParams[0]).get()
          chainList = ['Any'] + inputObj.getChainNames()
        except Exception as e:
          print("ERROR: ", e)
          return

        finalChainList = []
        for i in chainList:
          finalChainList.append(String('Chain: ' + i))
        provider = ListTreeProviderString(finalChainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                              "Select one of the chains (model, chain, "
                              "number of chain residues)")
        form.setVar(outputParam[0], dlg.values[0].get())

SelectChainWizardGromacs().addTarget(protocol=GromacsMDSimulation,
                                     targets=['restrainChain'],
                                     inputs=['gromacsSystem'],
                                     outputs=['restrainChain'])

class SelectResidueWizardGromacs(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getPrevResidueRest(self, groups):
        prev = []
        for idx, val in groups.items():
            if 'ResidueRestraint' in val:
                prev.append(val)
        return prev

    def parseTopoResidues(self, topFile, chains, chainLengths):
        '''Returns the residues as {chainName: {resIdx: resType, ...}, ...}'''
        resDic = {c: {} for c in chains}
        chainIdx = 0
        current_chain = chains[chainIdx]
        inAtoms = False

        with open(topFile) as f:
            for line in f:
                line = line.strip()

                # Wait until we hit the [ atoms ] section
                if not inAtoms:
                    if line.startswith('[ atoms ]'):
                        inAtoms = True
                    continue

                # Stop parsing if we hit an empty line after [ atoms ]
                if not line:
                    break

                # Parse residue lines
                if line.startswith('; residue'):
                    parts = line.split()
                    resIdx = int(parts[2])
                    resType = parts[3]

                    # Assign residue to the current chain
                    resDic[current_chain][resIdx] = resType

                    # If the current chain has all its expected residues, move to the next chain
                    if len(resDic[current_chain]) == chainLengths[chainIdx]:
                        chainIdx += 1
                        # Update current_chain if there are still chains left to process
                        if chainIdx < len(chains):
                            current_chain = chains[chainIdx]

        return resDic

    def getResidues(self, inputObj):
        '''Returns the residues in the selected chain as [[resIdx, resType], ...., [resIdx, resType]]'''
        chains = inputObj.getChainNames()
        chainLengths = inputObj.getChainLengths()
        topFile = inputObj.getTopologyFile()

        resDic = self.parseTopoResidues(topFile, chains, chainLengths)
        return resDic

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        inputObj = getattr(protocol, inputParams[0]).get()
        chainStr = getattr(protocol, inputParams[1]).get().split()[1]

        resDic = self.getResidues(inputObj)

        ALL_RES = 'All'
        residueList = [String(ALL_RES)]
        for chain in resDic:
            if chainStr in ['Any', chain]:
                for resIdx in resDic[chain]:
                    resType = resDic[chain][resIdx]
                    residueList.append(String('{"chain": "%s", "index": %d, "residue": "%s"}' %
                                              (chain, resIdx, resType)))

        provider = ListTreeProviderString(residueList)
        dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                              "Select one or several residues (residue number, "
                              "residue name)")

        groups = gromacsPlugin.parseIndexFile(protocol, gromacsPlugin.ensureIndexFile(protocol))
        restraintStr = 'ResidueRestraint_{}: '.format(len(self.getPrevResidueRest(groups)))

        if len(dlg.values) == 1 and dlg.values[0].get() == ALL_RES:
            for chain in resDic:
                if chainStr in ['Any', chain]:
                    restraintStr += '{}_{} '.format(chain, ALL_RES)

        else:
            idxsDic = {ch: [] for ch in inputObj.getChainNames()}
            for val in dlg.values:
                if val.get() != ALL_RES:
                    idxsDic[json.loads(val.get())['chain']].append(json.loads(val.get())['index'])

            for ch in idxsDic:
                if idxsDic[ch]:
                    idxGroups = groupConsecutiveIdxs(idxsDic[ch])
                    idxStr = ['{}-{}'.format(ig[0], ig[-1]) if len(ig) > 1 else str(ig[0]) for ig in idxGroups]
                    restraintStr += '{}_{} '.format(ch, ','.join(idxStr))

        form.setVar(outputParam[0], restraintStr)

SelectResidueWizardGromacs().addTarget(protocol=GromacsMDSimulation,
                                       targets=['restrainResidue'],
                                       inputs=['gromacsSystem', 'restrainChain'],
                                       outputs=['restrainResidue'])



class AddROIRestraintWizard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getSelectedROI(self, roiPointers, roiIdx, roiName):
        roi = None
        for roi in roiPointers[roiIdx].get():
            if roi.__str__() == roiName:
                break
        if not roi:
            print('Selected ROI not found')
        return roi

    def createROIRestraintIndex(self, system, roi, protocol):
        inIndex, outIndex = gromacsPlugin.ensureIndexFile(protocol), gromacsPlugin.getCustomIndexFile(protocol)

        #  Creating index for atoms in ROI
        atomGroups = groupConsecutiveIdxs(roi.getDecodedCAtoms())
        inCommand = ' | '.join(['a {}-{}'.format(ag[0], ag[-1]) for ag in atomGroups])
        gromacsPlugin.createIndexFile(protocol, system, inputCommands=[inCommand], inIndex=inIndex, outIndex=outIndex)
        groups = gromacsPlugin.parseIndexFile(protocol, outIndex)

        # Renaming ROI index to ROI name
        roiIdx, roiName = list(groups.keys())[-1], roi.__str__().replace(' ', '_')
        inCommand = 'name {} {}'.format(roiIdx, roiName)
        gromacsPlugin.createIndexFile(protocol, system, inputCommands=[inCommand], inIndex=outIndex, outIndex=outIndex)
        return roiIdx, roiName

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        prevPointers = getattr(protocol, outputParam[0])
        prevIds = protocol.getPrevPointersIds(prevPointers)
        newSet = getattr(protocol, inputParams[1]).get()
        newId = newSet.getObjId()

        if newId not in prevIds:
            newIndex = len(prevPointers)
            prevPointers.append(Pointer(newSet))
        else:
            newIndex = prevIds.index(newId)
        form.setVar(outputParam[0], prevPointers)

        roiName = getattr(protocol, inputParams[2]).get()
        roi = self.getSelectedROI(prevPointers, newIndex, roiName)
        system = getattr(protocol, inputParams[0]).get()
        _, roiName = self.createROIRestraintIndex(system, roi, protocol)

AddROIRestraintWizard().addTarget(protocol=GromacsMDSimulation,
                                  targets=['restraintROIInfo'],
                                  inputs=['gromacsSystem', 'restrainROIs', 'restrainROI'],
                                  outputs=['inputPointers'])

class AddResidueRestraintWizard(SelectResidueWizardGromacs):
    _targets, _inputs, _outputs = [], {}, {}

    def parseTopoAtoms(self, topFile, chains, chainLengths):
        '''Returns the atoms grouped by residue as
        {chainName: {resIdx: [atomId1, atomId2, ...]}, ...}'''
        resDic = {c: {} for c in chains}
        chain_limits = iter(zip(chains, chainLengths))
        curr_chain, target_len = next(chain_limits, (None, 0))
        curr_res = None

        with open(topFile) as f:
            for line in f:
                if line.startswith('[ atoms ]'):
                    break

            # Parse the atoms block
            for line in f:
                if not line.strip():
                    break  # End of the block

                if line.startswith('; residue') and curr_chain:
                    curr_res = int(line.split()[2])
                    resDic[curr_chain][curr_res] = []

                    # Move to the next chain only when this chain reaches its expected residue count
                    if len(resDic[curr_chain]) == target_len:
                        curr_chain, target_len = next(chain_limits, (None, 0))

                elif line.strip() and not line.startswith(';') and curr_chain and curr_res is not None:
                    atom_id = int(line.split()[0])
                    resDic[curr_chain][curr_res].append(atom_id)
        return resDic

    def getResDic(self, resStr):
        resDic = {}
        for chainRes in resStr.split()[1:]:
            print('ChainRes: ', chainRes)
            ch, res = chainRes.split('_')
            resDic[ch] = res
        return resDic

    def expandInterval(self, interval):
        if '-' in interval:
            first, last = interval.split('-')
            inds = list(range(int(first), int(last) + 1))
        else:
            inds = [int(interval)]
        return inds

    def createResiduesRestraintIndex(self, system, resDic, protocol):
        chains = system.getChainNames()
        chainLengths = system.getChainLengths()
        topFile = system.getTopologyFile()
        atomsDic = self.parseTopoAtoms(topFile, chains, chainLengths)

        inIndex, outIndex = gromacsPlugin.ensureIndexFile(protocol), gromacsPlugin.getCustomIndexFile(protocol)

        # Mapping residues to atoms index
        atomList = []
        for ch in resDic:
            if resDic[ch] == 'All':
                for resIdx in atomsDic[ch]:
                    atomList += atomsDic[ch][resIdx]

            elif resDic[ch]:
                for resGroup in resDic[ch].split(','):
                    for resIdx in self.expandInterval(resGroup):
                        atomList += atomsDic[ch][resIdx]

        atomGroups = groupConsecutiveIdxs(atomList)
        inCommand = ' | '.join(['a {}-{}'.format(ag[0], ag[-1]) for ag in atomGroups])

        gromacsPlugin.createIndexFile(protocol, system, inputCommands=[inCommand], inIndex=inIndex, outIndex=outIndex)
        groups = gromacsPlugin.parseIndexFile(protocol, outIndex)

        # Renaming group
        restIdx, restName = list(groups.keys())[-1], 'ResidueRestraint_{}'.format(len(self.getPrevResidueRest(groups)))
        inCommand = 'name {} {}'.format(restIdx, restName)
        gromacsPlugin.createIndexFile(protocol, system, inputCommands=[inCommand], inIndex=outIndex, outIndex=outIndex)
        return restIdx, restName

    def show(self, form, *params):
        inputParams, _ = self.getInputOutput(form)
        protocol = form.protocol

        resStr = getattr(protocol, inputParams[1]).get()
        system = getattr(protocol, inputParams[0]).get()
        resDic = self.getResDic(resStr)
        print(resDic)
        self.createResiduesRestraintIndex(system, resDic, protocol)

AddResidueRestraintWizard().addTarget(protocol=GromacsMDSimulation,
                                  targets=['restraintResidueInfo'],
                                  inputs=['gromacsSystem', 'restrainResidue'],
                                  outputs=[])

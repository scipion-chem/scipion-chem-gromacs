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
import json
from pyworkflow.gui import ListTreeProviderString, dialog

from pwem.objects import Pointer, String

from pwchem.wizards import AddElementSummaryWizard, DeleteElementWizard, VariableWizard, SelectElementWizard, \
    WatchElementWizard
from pwchem.utils import groupConsecutiveIdxs

from ..protocols.protocol_MD_simulation import *

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
        inputParam, outputParam = self.getInputOutput(form)

        system = getattr(protocol, inputParam[0]).get()

        outIndex = protocol.getCustomIndexFile()
        if os.path.exists(outIndex):
            groups = protocol.parseIndexFile(outIndex)
        else:
            groups = protocol.createIndexFile(system, inIndex=None, outIndex=outIndex)

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
        inputParam, outputParam = self.getInputOutput(form)

        system = getattr(protocol, inputParam[0]).get()

        inIndex, outIndex = None, protocol.getCustomIndexFile()
        if os.path.exists(outIndex):
            inIndex = outIndex

        inCommand = getattr(protocol, inputParam[1]).get()
        inCommand = ' '.join(protocol.translateNamesToIndexGroup(inCommand.split()))

        groups = protocol.createIndexFile(system, inputCommands=[inCommand, ''], inIndex=inIndex, outIndex=outIndex)


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
          finalChainList.append(pwobj.String('Chain: ' + i))
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

    def parseTopoResidues(self, topFile, chains):
        '''Returns the residues as {chainName: {resIdx: resType, resIdx2: resType2}, ....,
        chainName2: {}, ...}'''
        inAtoms, prevIdx = False, None
        resDic, chainIdx = {c: {} for c in chains}, 0
        with open(topFile) as f:
            for line in f:
                chainStr = chains[chainIdx]
                if inAtoms:
                    if line.startswith('; residue'):
                        resIdx, resType = line.split()[2:4]
                        if not prevIdx or int(resIdx) - 1 == prevIdx:
                            resDic[chainStr][int(resIdx)] = resType
                        else:
                            chainIdx += 1
                            resDic[chains[chainIdx]][int(resIdx)] = resType
                        prevIdx = int(resIdx)

                if not inAtoms and line.startswith('[ atoms ]'):
                    inAtoms = True
                elif inAtoms and not line.strip():
                    break
        return resDic

    def getResidues(self, inputObj):
        '''Returns the residues in the selected chain as [[resIdx, resType], ...., [resIdx, resType]]'''
        chains = inputObj.getChainNames()
        topFile = inputObj.getTopologyFile()

        resDic = self.parseTopoResidues(topFile, chains)
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

        groups = protocol.parseIndexFile(protocol.getCustomIndexFile())
        restraintStr = 'ResidueRestraint_{}: '.format(len(self.getPrevResidueRest(groups)))

        if len(dlg.values) == 1 and dlg.values[0].get() == ALL_RES:
            for chain in resDic:
                if chainStr in ['Any', chain]:
                    restraintStr += '{}_{} '.format(chain, ALL_RES)

        else:
            idxsDic = {ch: [] for ch in inputObj.getChainNames()}
            for val in dlg.values:
                if not val.get() == ALL_RES:
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
        inIndex, outIndex = None, protocol.getCustomIndexFile()
        if os.path.exists(outIndex):
            inIndex = outIndex

        #  Creating index for atoms in ROI
        atomGroups = groupConsecutiveIdxs(roi.getDecodedCAtoms())
        inCommand = ' | '.join(['a {}-{}'.format(ag[0], ag[-1]) for ag in atomGroups])
        groups = protocol.createIndexFile(system, inputCommands=[inCommand, ''], inIndex=inIndex, outIndex=outIndex)

        # Renaming ROI index to ROI name
        roiIdx, roiName = list(groups.keys())[-1], roi.__str__().replace(' ', '_')
        inCommand = 'name {} {}'.format(roiIdx, roiName)
        protocol.createIndexFile(system, inputCommands=[inCommand, ''], inIndex=outIndex, outIndex=outIndex)
        return roiIdx, roiName

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        prevPointers = getattr(protocol, outputParam[0])
        prevIds = protocol.getPrevPointersIds(prevPointers)
        newSet = getattr(protocol, inputParams[1]).get()
        newId = newSet.getObjId()

        if not newId in prevIds:
            newIndex = len(prevPointers)
            prevPointers.append(Pointer(newSet))
        else:
            newIndex = prevIds.index(newId)
        form.setVar(outputParam[0], prevPointers)

        roiName = getattr(protocol, inputParams[2]).get()
        roi = self.getSelectedROI(prevPointers, newIndex, roiName)
        system = getattr(protocol, inputParams[0]).get()
        roiIdx, roiName = self.createROIRestraintIndex(system, roi, protocol)

AddROIRestraintWizard().addTarget(protocol=GromacsMDSimulation,
                                  targets=['restraintROIInfo'],
                                  inputs=['gromacsSystem', 'restrainROIs', 'restrainROI'],
                                  outputs=['inputPointers'])

class AddResidueRestraintWizard(SelectResidueWizardGromacs):
    _targets, _inputs, _outputs = [], {}, {}

    def parseTopoAtoms(self, topFile, chains):
        '''Returns the residues as
        {chainName: {resIdx: atoms, resIdx2: atoms2}, chainName2: {resIdx: atoms, ...}, ...}'''
        inAtoms, prevIdx = False, None
        resDic, chainIdx = {c: {} for c in chains}, 0
        with open(topFile) as f:
            for line in f:
                chainStr = chains[chainIdx]
                if inAtoms:
                    if line.startswith('; residue'):
                        resIdx = line.split()[2]
                        if not prevIdx or int(resIdx) - 1 == prevIdx:
                            resDic[chainStr][int(resIdx)] = []
                        else:
                            chainIdx += 1
                            resDic[chains[chainIdx]][int(resIdx)] = []
                        prevIdx = int(resIdx)
                    elif line.strip() and not line.startswith(';'):
                        resDic[chainStr][prevIdx].append(int(line.split()[0]))

                if not inAtoms and line.startswith('[ atoms ]'):
                    inAtoms = True
                elif inAtoms and not line.strip():
                    break
        return resDic

    def getResDic(self, resStr):
        resDic = {}
        for chain_res in resStr.split()[1:]:
            print('Chain_res: ', chain_res)
            ch, res = chain_res.split('_')
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
        topFile = system.getTopologyFile()
        atomsDic = self.parseTopoAtoms(topFile, chains)

        inIndex, outIndex = None, protocol.getCustomIndexFile()
        if os.path.exists(outIndex):
            inIndex = outIndex

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
        groups = protocol.createIndexFile(system, inputCommands=[inCommand, ''], inIndex=inIndex, outIndex=outIndex)

        # Renaming group
        restIdx, restName = list(groups.keys())[-1], 'ResidueRestraint_{}'.format(len(self.getPrevResidueRest(groups)))
        inCommand = 'name {} {}'.format(restIdx, restName)
        groups = protocol.createIndexFile(system, inputCommands=[inCommand, ''], inIndex=outIndex, outIndex=outIndex)
        return restIdx, restName

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        resStr = getattr(protocol, inputParams[1]).get()
        system = getattr(protocol, inputParams[0]).get()
        resDic = self.getResDic(resStr)
        restIdx, restName = self.createResiduesRestraintIndex(system, resDic, protocol)

AddResidueRestraintWizard().addTarget(protocol=GromacsMDSimulation,
                                  targets=['restraintResidueInfo'],
                                  inputs=['gromacsSystem', 'restrainResidue'],
                                  outputs=[])


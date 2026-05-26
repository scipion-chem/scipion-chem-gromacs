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
import json, os, re
from Bio import PDB
from pyworkflow.gui import ListTreeProviderString, dialog
import pyworkflow.object as pwobj

from pwem.objects import Pointer, String

from pwchem.utils import groupConsecutiveIdxs
from pwchem.utils import pdbFromASFile
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
        resDic = {c: {} for c in chains}

        # Create the iterator and set the initial state
        chainLimits = iter(zip(chains, chainLengths))
        currChain, targetLen = next(chainLimits, (None, 0))

        inAtoms = False

        with open(topFile) as f:
            for line in f:
                line = line.strip()

                if not line:
                    continue

                # Only parse inside [ atoms ]
                if line.startswith('['):
                    inAtoms = 'atoms' in line
                    continue

                # Parse residues
                if inAtoms and line.startswith('; residue') and currChain:

                    # 1. Switch to the next chain if the current one is full
                    if len(resDic[currChain]) == targetLen:
                        currChain, targetLen = next(chainLimits, (None, 0))

                    # 2. Extract and assign the residue data
                    parts = line.split()
                    if len(parts) >= 4 and currChain:
                        resIdx = int(parts[2])
                        resType = parts[3]

                        resDic[currChain][resIdx] = resType
        return resDic

    def getResidues(self, inputObj):
        """Returns all residues as {chain: {resIdx: resType}}."""
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
        chainLimits = iter(zip(chains, chainLengths))
        currChain, targetLen = next(chainLimits, (None, 0))
        currRes = None

        with open(topFile) as f:
            for line in f:
                if line.startswith('[ atoms ]'):
                    break

            # Parse the atoms block
            for line in f:
                if not line.strip():
                    break  # End of the block

                if line.startswith('; residue') and currChain:
                    # Move to the next chain BEFORE processing this residue
                    # if we've already reached the expected residue count for current chain
                    if len(resDic[currChain]) == targetLen:
                        currChain, targetLen = next(chainLimits, (None, 0))

                    # Now assign this residue to the correct chain
                    currRes = int(line.split()[2])
                    if currChain:  # Make sure we still have a chain to assign to
                        resDic[currChain][currRes] = []

                elif line.strip() and not line.startswith(';') and currChain and currRes is not None:
                    atomId = int(line.split()[0])
                    resDic[currChain][currRes].append(atomId)

        return resDic

    def getResDic(self, resStr):
        resDic = {}
        for chainRes in resStr.split()[1:]:
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
        self.createResiduesRestraintIndex(system, resDic, protocol)

AddResidueRestraintWizard().addTarget(protocol=GromacsMDSimulation,
                                  targets=['restraintResidueInfo'],
                                  inputs=['gromacsSystem', 'restrainResidue'],
                                  outputs=[])


class SelectSSBondWIzard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParams, outputParams = self.getInputOutput(form)
        protocol = form.protocol

        # Detect SS bonds by running pdb2gmx in dry-run mode
        outPath = os.path.abspath(protocol.getProject().getTmpPath('SS_info.log'))
        print(getattr(protocol, inputParams[0]).get())

        inputStruct = os.path.abspath(getattr(protocol, inputParams[0]).get().getFileName())
        base, ext = os.path.splitext(os.path.basename(inputStruct))

        # Only convert if the input is not already a PDB
        if ext.lower() != '.pdb':
            tmpPdbPath = os.path.abspath(protocol.getProject().getTmpPath(f'{base}_temp.pdb'))
            pdbFile = pdbFromASFile(inputStruct, tmpPdbPath)
        else:
            pdbFile = inputStruct

        ssBonds = self.detectSSBonds(pdbFile, outPath)
        print(ssBonds)

        if not ssBonds:
            dialog.showInfo("Disulfide Bonds",
                            "No potential disulfide bonds detected in structure.",
                            form.root)
            return

        # 1. Wrap the bond labels into pyworkflow String objects
        bondStringList = [String(bond['label']) for bond in ssBonds]

        # 2. Create the provider and the ListDialog
        provider = ListTreeProviderString(bondStringList)
        dlg = dialog.ListDialog(form.root, "Select Disulfide Bonds", provider,
                                "Select which disulfide bonds to form\n(Ctrl+Click or Shift+Click for multiple):")

        # 3. Process the selections if the user clicked OK and made a choice
        if dlg.values:
            selectedIndices = []
            selectedLabels = [val.get() for val in dlg.values]

            # Map the selected labels back to their original 0-based indices
            for i, bond in enumerate(ssBonds):
                if bond['label'] in selectedLabels:
                    selectedIndices.append(str(i))

            if selectedIndices:
                form.setVar(outputParams[0], ','.join(selectedIndices))

    def detectSSBonds(self, inputStructure, outputFile, waterff='spc', mainff='amber03'):
        """
        Run pdb2gmx with -ss to count how many SS bonds GROMACS detects.
        Parse the output to extract the bond proposals.
        """
        maxBonds = self.maximumSSbonds(inputStructure)
        params = f'pdb2gmx -f {inputStructure} -o test.gro -water {waterff} -ff {mainff} -merge all -ss -ignh' \
               f' > {outputFile} 2>&1'

        printfValues =['n'] * maxBonds
        gromacsPlugin.runGromacsPrintfViewer(printfValues, params, cwd=os.path.dirname(outputFile))

        # Parse the output to extract SS bond proposals
        bonds = self.parseSSBondOutput(outputFile)
        return bonds

    def maximumSSbonds(self, inputStructure):
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure("protein", inputStructure)
        cysCount = sum(1 for residue in structure.get_residues()
                        if residue.get_resname() in ['CYS', 'CYX'])
        # Maximum possible bonds is = n*(n-1)/2
        max_bonds = (cysCount * (cysCount - 1)) // 2
        return max_bonds

    def parseSSBondOutput(self, outputFile):
        """
        Parse GROMACS output to extract SS bond proposals.
        Returns list of dicts: [{'cys1': 'CYS-3', 'cys2': 'CYS-40', 'label': 'CYS-3 and CYS-40']
        """
        with open(outputFile, 'r') as f:
            content = f.read()
        bonds = []
        # Parse bond proposals
        pattern = r'Link\s+(CYS-\d+)\s+SG-\d+\s+and\s+(CYS-\d+)\s+SG-\d+\s+\(y/n\)\s*\?'

        for match in re.finditer(pattern, content):
            cys1, cys2 = match.groups()

            # Extract residue numbers for distance lookup
            res1 = int(cys1.split('-')[1])
            res2 = int(cys2.split('-')[1])

            bonds.append({
                'cys1': cys1,
                'cys2': cys2,
                'label': f"{cys1} and {cys2}"
            })

        return bonds

SelectSSBondWIzard().addTarget(protocol=GromacsSystemPrep,
                                  targets=['selectSSBonds'],
                                  inputs=['inputStructure'],
                                  outputs=['selectSSBonds'])
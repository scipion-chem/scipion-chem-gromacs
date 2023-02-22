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
import pyworkflow.wizard as pwizard
from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String

from pwchem.wizards import AddElementSummaryWizard, DeleteElementWizard, VariableWizard

from ..protocols.protocol_MD_simulation import *

AddElementSummaryWizard().addTarget(protocol=GromacsMDSimulation,
                             targets=['insertStep'],
                             inputs=['insertStep'],
                             outputs=['workFlowSteps', 'summarySteps'])

DeleteElementWizard().addTarget(protocol=GromacsMDSimulation,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])

class GromacsWatchRelaxStepWizard(pwizard.Wizard):
    """Watch the parameters of the step of the workflow defined by the index"""
    _targets = [(GromacsMDSimulation, ['watchStep'])]

    def show(self, form, *params):
        protocol = form.protocol
        watchStep = protocol.watchStep.get().strip()
        try:
            index = int(watchStep)
            if protocol.countSteps() >= index > 0:
                workSteps = protocol.workFlowSteps.get().split('\n')
                msjDic = eval(workSteps[index - 1])
                for pName in msjDic:
                    if pName in protocol._paramNames:
                        form.setVar(pName, msjDic[pName])
                    elif pName in protocol._enumParamNames:
                        if pName == 'integrator':
                            idx = protocol._integrators.index(msjDic[pName])
                        elif pName == 'ensemType':
                            idx = protocol._ensemTypes.index(msjDic[pName])
                        elif pName == 'thermostat':
                            idx = protocol._thermostats.index(msjDic[pName])
                        elif pName == 'barostat':
                            idx = protocol._barostats.index(msjDic[pName])
                        elif pName == 'coupleStyle':
                            idx = protocol._coupleStyle.index(msjDic[pName])
                        elif pName == 'restraints':
                            idx = protocol._restraintTypes.index(msjDic[pName])
                        form.setVar(pName, idx)
        except:
            print('Index "{}" not recognized as integer for watch step'.format(watchStep))


class GromacsCheckIndexWizard(VariableWizard):
    """Watch the groups contained in the input Gromacs system"""
    _targets, _inputs, _outputs = [], {}, {}

    def createGroupsFile(self, system, inIndex=None, outIndex='/tmp/indexes.ndx', outFile='/tmp/indexGroups.txt',
                           inputCommands=['q']):
        outDir = os.path.dirname(outFile)
        inIndex = ' -n {}'.format(inIndex) if inIndex else ''
        command = 'make_ndx -f {}{} -o {} > {}'.format(system.getSystemFile(), inIndex, outIndex, outFile)

        if not inputCommands[-1] == 'q':
            inputCommands.append('q')
        gromacsPlugin.runGromacsPrintf(printfValues=inputCommands, args=command, cwd=outDir)
        return outFile

    def show(self, form, *params):
        protocol = form.protocol
        inputParam, outputParam = self.getInputOutput(form)

        system = getattr(protocol, inputParam[0]).get()
        outFile = protocol.getCustomGroupsFile()

        if not os.path.exists(outFile):
            inIndex, outIndex = None, protocol.getCustomIndexFile()
            if os.path.exists(outIndex):
                inIndex = outIndex

            outFile = self.createGroupsFile(system, inIndex=inIndex, outIndex=protocol.getCustomIndexFile(),
                                            outFile=outFile)

        groups = protocol.parseGroupsFile(outFile)

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

        outFile = protocol.getCustomGroupsFile()
        inIndex, outIndex = None, protocol.getCustomIndexFile()
        if os.path.exists(outIndex):
            inIndex = outIndex

        inCommand = getattr(protocol, inputParam[1]).get()
        inCommand = ' '.join(protocol.translateNamesToIndexGroup(inCommand.split()))

        groups = self.createGroupsFile(system, inputCommands=[inCommand], inIndex=inIndex, outIndex=outIndex,
                                         outFile=outFile)
        # Needs to be run twice. Seems like first with commands does not create the groups file, but it does create
        # the index file
        self.createGroupsFile(system, inIndex=inIndex, outFile=outFile)

        # form.setVar(outputParam[0], inCommand)


GromacsCustomIndexWizard().addTarget(protocol=GromacsMDSimulation,
                               targets=['restraintCommand'],
                               inputs=['gromacsSystem', 'restraintCommand'],
                               outputs=['restraints'])


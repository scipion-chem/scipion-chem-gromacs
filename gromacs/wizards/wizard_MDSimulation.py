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
from ..protocols.protocol_MD_simulation import *
import pyworkflow.wizard as pwizard
from pwchem.wizards import AddElementWizard, DeleteElementWizard

AddElementWizard().addTarget(protocol=GromacsMDSimulation,
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
                        elif pName == 'restrains':
                            idx = protocol._restrainTypes.index(msjDic[pName])
                        form.setVar(pName, idx)
        except:
            print('Index "{}" not recognized as integer for watch step'.format(watchStep))


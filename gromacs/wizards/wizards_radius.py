# **************************************************************************
# *
# * Authors:     Jose Gutierrez (jose.gutierrez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import numpy as np

from pwem.emlib.image import ImageHandler
from pwem.convert import AtomicStructHandler
from pwem.wizards.wizard import EmWizard

from pwem.wizards.wizards_3d.mask_structure_wizard import MaskStructureWizard
from gromacs.protocols.protocol_system_prep import GromacsSystemPrep


class ProtTestWizard_2(EmWizard):
    _targets = [(GromacsSystemPrep, ['radius_2', 'x', 'y', 'z'])]

    def show(self, form):
        protocol = form.protocol
        structure = protocol.UsePDBFile.get()
        if not structure:
            print('You must specify input structure')
            return
        plt = MaskStructureWizard(structure.getFileName())
        plt.initializePlot()
        form.setVar('radius_2', plt.radius)
        form.setVar('x', plt.origin[0])
        form.setVar('y', plt.origin[1])
        form.setVar('z', plt.origin[2])

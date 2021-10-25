# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Pedro Febrer Martinez (pedrofebrer98@gmail.com)
# *
# * your institution
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
from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pyworkflow.wizard import Wizard

"""
class GromacsPDB2GMX(Wizard):
    # Dictionary to target protocol parameters
    _targets = [(GromacsPDB2GMX, ['mainForceField']),
                (GromacsPDB2GMX, ['waterForceField'])]

    def show(self, form, *params):

        # This are the main force fields:
        PFF = [String("AMBER03"), String("AMBER94"),
                     String("AMBER96"), String("AMBER99"),
                     String("AMBER99SB"), String("AMBER99SB-ILDN"),
                     String("AMBERGS"), String("CHARMM27"),
                     String("GROMOS43a1"), String("GROMOS43a2"),
                     String("GROMOS45a3"), String("GROMOS53a5"),
                     String("GROMOS53a6"), String("GROMOS54a7"),
                     String("OPLS-AA")]


        # Get a data provider from the force fields to be used in the tree (dialog)
        provider = ListTreeProviderString(PFF)

        # Show the dialog
        dlg = dialog.ListDialog(form.root, "Available Main Force Fields", provider,
                                "Select one of the main force fields)")


        # Set the chosen value back to the form
        form.setVar('mainForceField', dlg.values[0].get())

    def show2(self, form, *params):
        WFF = [String("TIP3P"), String("TIP4P"),
                     String("TIP4PEW"), String("TIP5P"),
                     String("TIP5P (Ewald Sums)"), String("SPC"),
                     String("SPC/E"), String("None")]

        provider = ListTreeProviderString(WFF)

        dlg = dialog.ListDialog(form.root, "Available Water Force Fields", provider,
                                 "Select one of the water force fields")

        form.setVar('waterForceField', dlg.values[0].get())
"""
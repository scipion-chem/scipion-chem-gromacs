# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
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


import os
import glob

from tkinter import messagebox

from pyworkflow.protocol.params import LabelParam, EnumParam
from pyworkflow.viewer import DESKTOP_TKINTER, WEB_DJANGO, ProtocolViewer
from pwem.viewers import TableView, Chimera, ChimeraView
import pyworkflow.utils as pwutils

from gromacs.protocols.protocol_analysis import GromacsAnalysis
from tkinter import *

def errorWindow(tkParent, msg):
    try:
        # if tkRoot is null the error message may be behind
        # other windows
        messagebox.showerror("Error",  # bar title
                             msg,  # message
                             parent=tkParent)
    except:
        print(("Error:", msg))


class DARCViewer(ProtocolViewer):
    """ Viewer for Rosetta program DARC
    """
    _label = 'Gromacs IMG Viewer'
    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [GromacsAnalysis]
    DARC_SCORE = 'md_0_1_hbond.ndx'
    """
    if self.XtcFileSource == 1:
        new_file = basename((self.XtcFileLocal.get()).split(".")[0])
    else:
        new_file = basename(os.path.abspath(self.UseMDGroFile.get().getMD_xtc()).split(".")[0])
    RMS = "%s_rmsd.png" % (new_file)
    
    
    def __init__(self, **kwargs):
        ProtocolViewer.__init__(self, **kwargs)

        # Get the DARC score file
        self.DARC_SCORE = self.protocol._getPath(self.DARC_SCORE)

        # Get the complex files with the best ligand between all conformers given to darc
        self.darc_complex = glob.glob(os.path.join(self.protocol._getExtraPath(), "darc_complex", "*"))
        self.COMPLEX = [os.path.basename(i) for i in self.darc_complex]

        # Get the used and chosen ligands files between all conformers given to darc
        self.darc_ligands = glob.glob(os.path.join(self.protocol._getExtraPath(), "ligands", "*"))
        self.LIGANDS = [os.path.basename(i) for i in self.darc_ligands]
    
        if self.protocol.RmsAnalysis.get():
            self.RMS = glob.glob(
                os.path.join(self.protocol._getExtraPath(), )
            )
        
        if self.protocol.minimize_output.get():
            # Get the complex files with the best ligand between all conformers given to darc
            self.mini_complex = glob.glob(
                os.path.join(self.protocol._getExtraPath(), "minimization", "darc_complex", "*"))
            self.COMPLEX_mini = [os.path.basename(i) for i in self.mini_complex]
            # Get the used and chosen ligands files between all conformers given to darc
            self.mini_ligands = glob.glob(os.path.join(self.protocol._getExtraPath(), "minimization", "ligands", "*"))
            self.LIGANDS_mini = [os.path.basename(i) for i in self.mini_ligands]
    
        with open(self.DARC_SCORE) as score_file:
            results = score_file.readlines()
            results = [x.strip().split("\t") for x in results]
            self.headerList = results[0]
            self.dataList = results[1:]
    """
    def _defineParams(self, form):
        form.addSection(label='DARC results')
        # group = form.addGroup('Overall results')
        if self.protocol.RmsAnalysis.get():
            form.addParam('showfinalscore', LabelParam,
                          label="Table of ligand DARC scores",
                          help="Table of DARC Scores\n\n"
                               "TAG: \n"
                               "DARC Score:.\n"
                               "Total Energy: \n"
                               "Interface Energy: .\n"
                               "Interface HB: .\n"
                               "Total packstats: .\n"
                               "Interface_unsat: .\n"
                               "ThetaLig: .\n")
        else:
            form.addParam('showfinalscore', LabelParam,
                          label="Table of ligand DARC scores",
                          help="Table of DARC Scores\n\n"
                               "TAG: \n"
                               "DARC Score:.")
        """
        view = form.addGroup("Visualization in Chimerax")
        view.addParam('show_complex', LabelParam,
                      label="Visualize complex protein-ligand")

        view.addParam('complex_v', EnumParam,
                      choices=self.COMPLEX,
                      default=0,
                      label="Choose a complex: ",
                      help="")

        view.addParam('show_ligand', LabelParam,
                      label="Visualize only ligand conformation")

        view.addParam('ligand_v', EnumParam,
                      choices=self.LIGANDS,
                      default=0,
                      label="Choose a ligand: ",
                      help="")
                      
        if self.protocol.minimize_output.get():
            view_mini = form.addGroup("Visualizationof complexes after minimization")
            view_mini.addParam('show_complex_mini', LabelParam,
                               label="Visualize complex protein-ligand")

            view_mini.addParam('complex_v_mini', EnumParam,
                               choices=self.COMPLEX_mini,
                               default=0,
                               label="Choose a complex: ",
                               help="")

            view_mini.addParam('show_ligand_mini', LabelParam,
                               label="Visualize only ligand conformation")

            view_mini.addParam('ligand_v_mini', EnumParam,
                               choices=self.LIGANDS_mini,
                               default=0,
                               label="Choose a ligand: ",
                               help="")
        """
    def _getVisualizeDict(self): #Que cosa mostrar, la variable ligand_v_mini, ejecuta funcion visualize_v_mini
        return {
            #'complex_v': self._visualizeComplex,
            #'ligand_v': self._visualizeLigand,
            #'complex_v_mini': self._visualize_miniComplex,
            #'ligand_v_mini': self._visualize_miniLigand,
            'showfinalscore': self._visualizeFinalResults}

    def _visualizeFinalResults(self, e=None):
        """ Create a formatted table with the DARC results, sorted by DARC score"""
        root = Tk()
        img = PhotoImage(file='/home/pedro/ScipionUserData/projects/Hinge_Region_GSK3/Runs/001009_GromacsAnalysis/md_0_1_hydrogen2-bond-intra-protein.png')
        label = Label(root, image=img)
        root.mainloop()
        print("Hello World")

    def _visualizeComplex(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.COMPLEX[self.complex_v.get()]
        path = os.path.abspath(os.path.join(self.protocol._getExtraPath(), "darc_complex", basename))
        self.displayChimera(file_path=path)

    def _visualizeLigand(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.LIGANDS[self.ligand_v.get()]
        path = os.path.abspath(os.path.join(self.protocol._getExtraPath(), "ligands", basename))
        self.displayChimera(file_path=path)

    def _visualize_miniComplex(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.COMPLEX_mini[self.complex_v_mini.get()]
        path = os.path.abspath(os.path.join(self.protocol._getExtraPath(), "minimization", "darc_complex", basename))
        self.displayChimera(file_path=path)

    def _visualize_miniLigand(self, e=None):
        """Visualize a complex protein-small molecule in Chimera"""
        basename = self.LIGANDS_mini[self.ligand_v_mini.get()]
        path = os.path.abspath(os.path.join(self.protocol._getExtraPath(), "minimization", "ligands", basename))
        self.displayChimera(file_path=path)

    # ------------------- UTILS functions --------------------------

    def displayChimera(self, file_path=None):

        if file_path == None:
            errorWindow(self.getTkRoot(), "Any structure available")
            return

        pwutils.runJob(None, Chimera.getProgram(), file_path, env=Chimera.getEnviron())
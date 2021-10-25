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

import pwem.objects.data as data
import pyworkflow.object as pwobj

class TopolStruct(data.EMFile):
    """Represents a Topol.top file. """
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)

class PosreItp(data.EMFile):
    """Represents a posre.itp (Position restriction) file. """
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)

class GroFile(data.EMFile):
    """Represents a gro (Position restriction) file. """
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)

class GroFiles(data.EMObject):
    """Represents all gro files"""
    def __init__(self, posre_h = None, posre = None, topol= None, gro=None):
        super().__init__()
        self.Posre_h = pwobj.String(posre_h)
        self.Posre = pwobj.String(posre)
        self.Topol = pwobj.String(topol)
        self.Gro = pwobj.String(gro)
    def __str__(self):
        return ", ".join([self.getPosre_h(), self.getPosre(), self.getTopol(), self.getGro()])
    def getPosre_h(self):
        return self.Posre_h.get()
    def getPosre(self):
        return self.Posre.get()
    def getTopol(self):
        return self.Topol.get()
    def getGro(self):
        return self.Gro.get()

class EMgroFile(data.EMObject):
    """Represents an em from gromacs file."""
    def __init__(self, em_edr = None, em_log = None, em_gro = None, em_trr = None, em_tpr = None):
        super().__init__()
        self.EM_edr = pwobj.String(em_edr)
        self.EM_log = pwobj.String(em_log)
        self.EM_gro = pwobj.String(em_gro)
        self.EM_trr = pwobj.String(em_trr)
        self.EM_tpr = pwobj.String(em_tpr)
    def __str__(self):
        return ", ".join([self.getEM_log(), self.getEM_tpr(), self.getEM_trr(), self.getEM_gro(), self.getEM_edr()])
    def getEM_edr(self):
        return self.EM_edr.get()
    def getEM_log(self):
        return self.EM_log.get()
    def getEM_gro(self):
        return self.EM_gro.get()
    def getEM_trr(self):
        return self.EM_trr.get()
    def getEM_tpr(self):
        return self.EM_tpr.get()

class NVTgroFile(data.EMObject):
    """Represents an em from gromacs file."""
    def __init__(self, nvt_edr = None, nvt_log = None, nvt_gro = None, nvt_trr = None, nvt_tpr = None, nvt_cpt = None):
        super().__init__()
        self.NVT_cpt = pwobj.String(nvt_cpt)
        self.NVT_edr = pwobj.String(nvt_edr)
        self.NVT_log = pwobj.String(nvt_log)
        self.NVT_gro = pwobj.String(nvt_gro)
        self.NVT_trr = pwobj.String(nvt_trr)
        self.NVT_tpr = pwobj.String(nvt_tpr)
    def __str__(self):
        return ", ".join([self.getNVT_log(), self.getNVT_tpr(), self.getNVT_trr(), self.getNVT_gro(), self.getNVT_edr(), self.getNVT_cpt()])
    def getNVT_edr(self):
        return self.NVT_edr.get()
    def getNVT_log(self):
        return self.NVT_log.get()
    def getNVT_gro(self):
        return self.NVT_gro.get()
    def getNVT_trr(self):
        return self.NVT_trr.get()
    def getNVT_tpr(self):
        return self.NVT_tpr.get()
    def getNVT_cpt(self):
        return self.NVT_cpt.get()


class NPTgroFile(data.EMObject):
    """Represents an em from gromacs file."""
    def __init__(self, npt_edr = None, npt_log = None, npt_gro = None, npt_trr = None, npt_tpr = None, npt_cpt = None):
        super().__init__()
        self.NPT_cpt = pwobj.String(npt_cpt)
        self.NPT_edr = pwobj.String(npt_edr)
        self.NPT_log = pwobj.String(npt_log)
        self.NPT_gro = pwobj.String(npt_gro)
        self.NPT_trr = pwobj.String(npt_trr)
        self.NPT_tpr = pwobj.String(npt_tpr)
    def __str__(self):
        return ", ".join([self.getNPT_log(), self.getNPT_tpr(), self.getNPT_trr(), self.getNPT_gro(), self.getNPT_edr(), self.getNPT_cpt()])
    def getNPT_edr(self):
        return self.NPT_edr.get()
    def getNPT_log(self):
        return self.NPT_log.get()
    def getNPT_gro(self):
        return self.NPT_gro.get()
    def getNPT_trr(self):
        return self.NPT_trr.get()
    def getNPT_tpr(self):
        return self.NPT_tpr.get()
    def getNPT_cpt(self):
        return self.NPT_cpt.get()


class MDgroFile(data.EMObject):
    """Represents an em from gromacs file."""
    def __init__(self, md_edr = None, md_log = None, md_gro = None, md_xtc = None, md_tpr = None, md_cpt = None):
        super().__init__()
        self.MD_cpt = pwobj.String(md_cpt)
        self.MD_edr = pwobj.String(md_edr)
        self.MD_log = pwobj.String(md_log)
        self.MD_gro = pwobj.String(md_gro)
        self.MD_xtc = pwobj.String(md_xtc)
        self.MD_tpr = pwobj.String(md_tpr)
    def __str__(self):
        return ", ".join([self.getMD_log(), self.getMD_tpr(), self.getMD_xtc(), self.getMD_gro(), self.getMD_edr(), self.getMD_cpt()])
    def getMD_edr(self):
        return self.MD_edr.get()
    def getMD_log(self):
        return self.MD_log.get()
    def getMD_gro(self):
        return self.MD_gro.get()
    def getMD_xtc(self):
        return self.MD_xtc.get()
    def getMD_tpr(self):
        return self.MD_tpr.get()
    def getMD_cpt(self):
        return self.MD_cpt.get()

class AnalysisFile(data.EMObject):
    """Represents an em from gromacs file."""
    def __init__(self, trjconv = None, ndx = None, rms = None, rmsf = None):
        super().__init__()
        self.Trjconv = pwobj.String(trjconv)
        self.Ndx = pwobj.String(ndx)
        self.Rms = pwobj.String(rms)
        self.Rmsf = pwobj.String(rmsf)
    def __str__(self):
        return ", ".join([self.getTrjconv(), self.getNdx(), self.getRms(), self.getRmsf()])
    def getTrjconv(self):
        return self.Trjconv.get()
    def getNdx(self):
        return self.Ndx.get()
    def getRms(self):
        return self.Rms.get()
    def getRmsf(self):
        return self.Rmsf.get()

class TrjconvFile(data.EMObject):
    """Represents a Trjconv xtc file."""
    def __init__(self, trjconv = None):
        super().__init__()
        self.Trjconv = pwobj.String(trjconv)
    def __str__(self):
        return ", ".join([self.getTrjconv()])
    def getTrjconv(self):
        return self.Trjconv.get()

class XvgFile(data.EMObject):
    """Represents a Trjconv xtc file."""
    def __init__(self, xvg = None):
        super().__init__()
        self.Xvg = pwobj.String(xvg)
    def __str__(self):
        return ", ".join([self.getXvg()])
    def getXvg(self):
        return self.Xvg.get()


class SetEMgroFile(data.EMSet):
    ITEM_TYPE = EMgroFile
    FILE_TEMPLATE_NAME = 'setEMgroFile%s.sqlite'
    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)

class TopolFile(TopolStruct):
    def __init__(self, filename=None, **kwargs):
        TopolStruct.__init__(self, filename, **kwargs)
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

import pyworkflow.viewer as pwviewer
import pwem.viewers.views as views
from gromacs.objects import *
import pwem.viewers.showj as showj

class setIDView(views.ObjectView):
    def __init__(self, project, inputID, path, other="", viewParams={}):
        defaultViewParams={showj.MODE:"metadata"}
        defaultViewParams.update(viewParams)
        views.ObjectView.__init__(self, project, inputID, path, other, defaultViewParams)

class gromacsSetViewer(pwviewer.Viewer):
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [SetEMgroFile, TopolStruct]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views=[]
        cls=type(obj)
        if issubclass(cls, SetEMgroFile):
            views.append(setIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, TopolFile):
            views.append(self.textView([obj.getFileName()]))
        return views
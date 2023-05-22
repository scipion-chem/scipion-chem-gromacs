#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Daniel Del Hoyo Gómez (ddelhoyo@cnb.csic.es)
# # # *
# # # *
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import sys, os
import BioSimSpace as BSS


def parseParams(paramsFile):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split('::')
      paramsDic[key.strip()] = value.strip()
  return paramsDic


if __name__ == "__main__":
    pDic = parseParams(sys.argv[1])
    sysName = os.path.splitext(os.path.basename(pDic['groFile']))[0]

    system = BSS.IO.readMolecules([pDic['groFile'], pDic['topFile']])

    box, shell = None, None
    if pDic['sizeType'] == 'Image distance':
        boxSize = float(pDic['boxSize'])
        box, angles = BSS.Box.generateBoxParameters(pDic['boxType'], boxSize * BSS.Units.Length.nanometer)
    else:
        shell = float(pDic['shellDist']) * BSS.Units.Length.nanometer

    solvated = BSS.Solvent.solvate(pDic['wff'], molecule=system, box=box, shell=shell,
               is_neutral=bool(pDic['bssNeutral']), ion_conc=float(pDic['bssIonConc']))

    BSS.IO.saveMolecules("{}_solvated".format(sysName), solvated, ["gro87", "grotop"])


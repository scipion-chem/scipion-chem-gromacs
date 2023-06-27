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

EMIN, NVT, NPT = 'Energy min', 'NVT', 'NPT'
ensemDic = {EMIN: 'minimisation', NVT: 'equilibration', NPT: 'equilibration'}
saveFormats = ["gro87", "grotop"]

def parseParams(paramsFile):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split('::')
      paramsDic[key.strip()] = value.strip()
  return paramsDic


if __name__ == "__main__":
    pDic = parseParams(sys.argv[1])
    msjDic = eval(pDic['msjDic'])
    rIdxs = eval(pDic['rAtomIdxs'])
    stageDir = pDic['stageDir']

    sysName = os.path.splitext(os.path.basename(pDic['groFile']))[0]
    system = BSS.IO.readMolecules([pDic['groFile'], pDic['topFile']])

    if msjDic['ensemType'] == EMIN:
      protocol = BSS.Protocol.Minimisation(steps=msjDic['nStepsMin'],
                                           restraint=rIdxs, force_constant=msjDic['restraintForce'])

    elif msjDic['ensemType']  in [NVT, NPT]:
      if msjDic['saveTrj']:
        saveFormats.append('trr')

      tCoup = msjDic['timeNeigh'] if msjDic['tempCouple'] == -1 else msjDic['tempCouple']
      stepsSave = int(msjDic['trajInterval'] / msjDic['timeStep'])


      protocol = BSS.Protocol.Equilibration(timestep=msjDic['timeStep'] * BSS.Units.Time.picosecond,
                                            runtime=msjDic['simTime'] * BSS.Units.Time.picosecond,
                                            temperature=msjDic['temperature'] * BSS.Units.Temperature.kelvin,
                                            thermostat_time_constant=tCoup * BSS.Units.Time.picosecond,
                                            restart_interval=stepsSave,
                                            restraint=rIdxs, force_constant=msjDic['restraintForce'])

    # Create a process object to run the simulation with AMBER.
    process = BSS.Process.Gromacs(system, protocol, work_dir=os.path.join(stageDir, 'workDir'))
    process.start()
    process.wait()

    # Get the minimised molecular system.
    minimised = process.getSystem()
    traj = process.getTrajectory()

    BSS.IO.saveMolecules('{}_min'.format(sysName), minimised, saveFormats)


# **************************************************************************
# *
# * Authors: Daniel Del Hoyo Gomez
# *          James M. Krieger
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

########################* Versions and Package dictionaries ###########################
V2020 = '2020.6'
V2021 = '2021.5'
V2022 = '2022.5'
V2023 = '2023.5'
V2024 = '2024.3'
GROMACS_VERSIONS = [V2021]
CMAKE_MINIMUM_VERSION = '3.16'

GROMACS_DIC = {'name': 'gromacs', 'version': V2021, 'home': 'GROMACS_HOME'}

MASTER = 'master'
V292 = '2.9.2'
V210B = '2.10b'
V210A = '2.10a'
PLUMED_VERSIONS = [MASTER, V292, V210B, V210A]
PLUMED_DIC = {'name': 'plumed', 'version': V210A, 'home': 'PLUMED_HOME'}

V200 = '2.0.0'
LIBTORCH_VERSIONS = [V200]
LIBTORCH_DIC = {'name': 'libtorch', 'version': V200, 'home': 'LIBTORCH_HOME'}

V01 = '0.1'
EMMIVOX_VERSIONS = [V01]
EMMIVOX_DIC = {'name': 'emmivox', 'version': V01, 'home': 'EMMIVOX_HOME'}

GROMACS_ENV_ACT = 'GROMACS_ENV_ACT'

"""
plumed 2.9.2 patching options:
1) gromacs-2020.7     3) gromacs-2022.5    5) gromacs-2024.2    7) namd-2.13         9) qespresso-5.0.2  11) qespresso-7.0
2) gromacs-2021.7     4) gromacs-2023.5    6) namd-2.12         8) namd-2.14        10) qespresso-6.2    12) qespresso-7.2

plumed 2.10b and master patching options:
 1) gromacs-2022.5
 2) gromacs-2023.5
 3) gromacs-2024.2
 4) namd-2.12
 5) namd-2.13
 6) namd-2.14
 7) qespresso-5.0.2
 8) qespresso-6.2
 9) qespresso-7.0
10) qespresso-7.2

plumed 2.10a patching options:
1) gromacs-2021.7     4) namd-2.12         7) qespresso-5.0.2  10) qespresso-7.2
2) gromacs-2022.5     5) namd-2.13         8) qespresso-6.2
3) gromacs-2023.2     6) namd-2.14         9) qespresso-7.0
"""

SCIPION_SOFTWARE = None

try:
    from pyworkflow import Config
    vars = Config.getVars()
    SCIPION_SOFTWARE = vars.get('SCIPION_SOFTWARE')
except Exception:
    pass

PATCH_DIC = {V292:  {V2020: 1, V2021: 2,
                     V2022: 3, V2023: 4,
                     V2024: 5},
            V210B:  {V2022: 1, V2023: 2,
                     V2024: 3},
            MASTER: {V2024: 3},
            V210A:  {V2021: 1, 
                     V2022: 2,
                     V2023: 3}}

########################################* ION NAMES #####################################

BR, CA, CL, CS, CU, CU2, F, I, K, LI, MG, NA, RB, ZN = 'BR-', 'CA2+', 'CL-', 'CS+', 'CU+', 'CU2+', 'F-', 'I-', 'K+', \
                                                       'LI+', 'MG2+', 'NA+', 'RB+', 'ZN2+'
ION_NAMES = ['BR', 'CA', 'CL', 'CS', 'CU', 'CU2', 'F', 'I', 'K', 'LI', 'MG', 'NA', 'RB', 'ZN']

###################################### MDP GENERATION #####################################
RESTR_STR = '''define = -DPOSRES_{}'''

TSTEP_EM = '''emtol = {}        ; Stop minimization when the maximum force < x kJ/mol/nm
emstep = {}          ; Minimization step size'''

TSTEP_EQ = '''dt = {}     ; Time Steps Size (ps)'''

DISP_CORR = '''DispCorr        = EnerPres  ; account for cut-off vdW scheme'''

OUTPUT_CONTROL = '''nstxout                 = {}       ; save coordinates every x * tStep ps
nstvout                 = {}       ; save velocities every x * tStep ps
nstlog                  = {}       ; update log file every x * tStep ps'''

BONDED_PARAMS = '''continuation            = {}       ; Restarting after another simulation 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy'''

ELECTROSTATICS = '''pme_order               = 4         ; cubic interpolation
coulombtype     = PME       ; Treatment of long range electrostatic interactions (PME: Particle Mesh Ewald)
fourierspacing          = 0.16      ; grid spacing for FFT'''

TEMP_SETTING = '''tcoupl                  = {}             ; thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
ref_t                   = {}     {}           ; reference temperature, one for each group, in K
tau_t                   = {}     {}           ; time constant, in ps
nsttcouple              = {}                  ; Frequency for temperature coupling
nstcomm                 = {}                  ; number of steps for center of mass motion removal'''

PRES_SETTING =  '''pcoupl                  = {}     ; Pressure coupling 
pcoupltype              = {}             ; uniform scaling of box vectors
ref_p                   = {}                  ; reference pressure, in bar
tau_p                   = {}                   ; time constant, in ps
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
nstpcouple              = {}                  ; Frequency for pressure coupling'''

VEL_GEN = '''gen_temp                = {}       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed'''


MDP_STR = '''; Parameters describing what to do, when to stop and what to save
; position restrains
{}  

; Run parameters
integrator  = {}         ; Algorithm 
nsteps      = {}         ; Maximum number of steps to perform
; Time Step Control
{}                       

; Nonbonded settings 
nstlist         = {}         ; Frequency to update the neighbor list and long range forces (fs)
cutoff-scheme   = Verlet    ; Buffered neighbor searching
rcoulomb        = 1.0       ; Short-range electrostatic cut-off (nm)
rvdw            = 1.0       ; Short-range Van der Waals cut-off (nm)
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions

; Disp Corr
{}

; Trajectory output control
{}

; Bond parameters
{}

; Electrostatics
{}

; Temperature coupling 
{}

; Pressure coupling 
{}

; Velocity generation
gen_vel                 = {}       ; assign velocities from Maxwell distribution
{}
'''
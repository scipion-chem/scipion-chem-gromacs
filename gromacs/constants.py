# **************************************************************************
# *
# * Authors: Daniel Del Hoyo Gomez
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
nstenergy               = {}       ; save energies every x * tStep ps
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


TCL_MD_STR = '''
mol addrep 0
mol new {%s} type {%s} first 0 last -1 step 1 waitfor 1
mol addfile {%s} type {%s} first 0 last -1 step 1 waitfor 1 0

mol color Name
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol selection protein
mol material Opaque
mol modrep 0 0

mol addrep 0
mol color Name
mol representation Points 1.000000
mol selection hetero within 3 of protein
mol material Opaque
mol modrep 1 0
'''

PML_MD_STR = '''load {}
load_traj {}
hide everything, not br. all within 3 of (byres polymer & name CA)
set movie_fps, 15
'''











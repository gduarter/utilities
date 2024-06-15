#!/bin/bash

# Author: Guilherme Duarte R. Matos
# Date: September 26, 2022

cat << EOF > 01.minimization.X.mdp
; Parameters describing what to do, when to stop and what to save
integrator      = steep     ; Algorithm (steep = steepest descent minimization)
emtol           = 1000.0    ; Stop minimization when the maximum force < 10.0 kJ/mol
emstep          = 0.01      ; Energy step size
nsteps          = 50000     ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1             ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type         = grid          ; Method to determine neighbor list (simple, grid)
rlist           = 1.2           ; Cut-off for making neighbor list (short range forces)
coulombtype     = PME           ; Treatment of long range electrostatic interactions
rcoulomb        = 1.2           ; long range electrostatic cut-off
vdwtype         = cutoff
vdw-modifier    = force-switch
rvdw-switch     = 1.0
rvdw            = 1.2           ; long range Van der Waals cut-off
pbc             = xyz           ; Periodic Boundary Conditions
DispCorr        = AllEnerPres
EOF

cat << EOF > 02.equil_nvt.X.mdp
title                   = Protein-ligand complex NVT equilibration 
define                  = -DPOSRES -DPOSRES_LIG ; position restrain the protein and ligand
; Run parameters
integrator              = sd        ; Langevin integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 0   
nstlog                  = 0   
nstxout-compressed      = 0   
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = hbonds    ; bonds to H are constrained 
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 20        ; largely irrelevant with Verlet
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb                = 1.2       ; short-range electrostatic cutoff (in nm)
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling
tcoupl                  = no        ; Langevin dynamics does not require coupling
tc-grps                 = System    ; two coupling groups - more accurate
tau_t                   = 2.0       ; time constant, in ps
ref_t                   = 298.15    ; reference temperature, one for each group, in K
; Pressure coupling
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = AllEnerPres
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 298.15    ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
ld_seed                 = -1
EOF

cat << EOF > 03.equil_npt.X.mdp
title                   = Protein-ligand complex NPT equilibration 
define                  = -DPOSRES -DPOSRES_LIG 
; Run parameters
integrator              = sd        ; Integrador de Langevin
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 0      
nstlog                  = 0      
nstxout-compressed      = 0      
; Bond parameters
continuation            = yes       ; continuação do NVT 
constraint_algorithm    = lincs      
constraints             = hbonds     
lincs_iter              = 1         
lincs_order             = 4         
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid     
nstlist                 = 20        
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2      
; Electrostatics
coulombtype             = PME       
rcoulomb                = 1.2
pme_order               = 4        
fourierspacing          = 0.16      
; Temperature coupling
tcoupl                  = no        
tc-grps                 = System    
tau_t                   = 2.0       
ref_t                   = 298.15    
; Pressure coupling
pcoupl                  = Berendsen    ; algoritmo do barostato para NPT
pcoupltype              = isotropic    ; caixa de simulacao varia uniformemente
tau_p                   = 2.0          ; constante de tempo
ref_p                   = 1.0          ; pressao de referencia em bar
compressibility         = 4.5e-5       ; compressibilidade isotermica da agua, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = AllEnerPres 
; Velocity generation
gen_vel                 = no        ; usar velocidades geradas da etapa NVT 
ld_seed                 = -1
EOF

cat << EOF > 04.equil_npt2.X.mdp
define                  = -DPOSRES -DPOSRES_LIG
; Run parameters
integrator              = sd        ; Integrador de Langevin
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstenergy               = 0
nstlog                  = 0
nstxout-compressed      = 0
; Bond parameters
continuation            = yes       ; continuação do NPT
constraint_algorithm    = lincs
constraints             = hbonds
lincs_iter              = 1
lincs_order             = 4
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
pme_order               = 4
fourierspacing          = 0.16
; Temperature coupling
tcoupl                  = no
tc-grps                 = System
tau_t                   = 2.0
ref_t                   = 298.15
; Pressure coupling
pcoupl                  = Parrinello-Rahman    ; algoritmo do barostato para NPT
pcoupltype              = isotropic            ; caixa de simulacao varia uniformemente
tau_p                   = 5.0                  ; constante de tempo
ref_p                   = 1.0                  ; pressao de referencia em bar
compressibility         = 4.5e-5               ; compressibilidade isotermica da agua, bar^-1
refcoord_scaling        = com                  ; permite o uso de restricoes de posicao com o barostato
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = AllEnerPres
; Velocity generation
gen_vel                 = no        ; usar velocidades geradas da etapa NPT
EOF

cat << EOF > 05.prod.X.mdp
; Run parameters
integrator              = sd          ; Integrador de Langevin
nsteps                  = 5000000     ; 10 ns
dt                      = 0.002       ; 2 fs
; Output control
nstenergy               = 1000
nstlog                  = 1000
nstxout-compressed      = 1000
; Bond parameters
continuation            = yes       ; continuação do NPT2
constraint_algorithm    = lincs
constraints             = hbonds
lincs_iter              = 1
lincs_order             = 4
; Neighbor searching and vdW
cutoff-scheme           = Verlet
ns_type                 = grid
nstlist                 = 20
rlist                   = 1.2
vdwtype                 = cutoff
vdw-modifier            = force-switch
rvdw-switch             = 1.0
rvdw                    = 1.2
; Electrostatics
coulombtype             = PME
rcoulomb                = 1.2
pme_order               = 4
fourierspacing          = 0.16
; Temperature coupling
tcoupl                  = no
tc-grps                 = System
tau_t                   = 2.0
ref_t                   = 298.15
; Pressure coupling
pcoupl                  = Parrinello-Rahman    ; algoritmo do barostato para NPT
pcoupltype              = isotropic            ; caixa de simulacao varia uniformemente
tau_p                   = 5.0                  ; constante de tempo
ref_p                   = 1.0                  ; pressao de referencia em bar
compressibility         = 4.5e-5               ; compressibilidade isotermica da agua, bar^-1
refcoord_scaling        = com                  ; permite o uso de restricoes de posicao com o barostato
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction is not used for proteins with the C36 additive FF
DispCorr                = AllEnerPres
; Velocity generation
gen_vel                 = no        ; usar velocidades geradas da etapa NPT
EOF



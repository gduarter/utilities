#!/bin/bash

###### Help Function ########
helpFunction(){
    echo -e "\tUsage: $0 -p pdb_code"
    exit 1
}

# Assign typed arguments to variables
while getopts "p:" opt
do
    case $opt in
        p ) PDB="$OPTARG";;
        ? ) helpFunction ;;
    esac
done

# Prints helpFuntion in case the number of parameters do not match what
# the script requires
if [ -z "${PDB}" ] 
then
    echo "You are misusing this script"
    helpFunction
fi

# define paths
root=$(pwd)
scriptdir=${root}/zzz.scripts
griddir=${root}/


cat << EOF > ${PDB}.multigrid.in
compute_grids                  yes
grid_spacing                   0.4
output_molecule                yes
contact_score                  no
chemical_score                 no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             9
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  ${prepdir}/${PDB}.rec.withH.mol2
box_file                       ${griddir}/${PDB}.box.pdb 
vdw_definition_file            ${DOCKHOME}/parameters/vdw_AMBER_parm99.defn
chemical_definition_file       ${DOCKHOME}/parameters/chem.defn
score_grid_prefix              ${PDB}.rec
receptor_out_file              ${PDB}.rec.grid.mol2
EOF


cat << EOF > ${PDB}.reference_multigrid.in
conformer_search_type                                        rigid
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             ${root}/${prepdir}/2nnq_lig_withH.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           no
grid_score_secondary                                         no
multigrid_score_primary                                      yes
multigrid_score_secondary                                    no
multigrid_score_rep_rad_scale                                1.0
multigrid_score_vdw_scale                                    1.0
multigrid_score_es_scale                                     1.0
multigrid_score_number_of_grids                              26                  
multigrid_score_grid_prefix0                                 ${PDB}.resid_  
multigrid_score_grid_prefix1                                 ${PDB}.resid_
multigrid_score_grid_prefix2                                 ${PDB}.resid_
multigrid_score_grid_prefix3                                 ${PDB}.resid_
multigrid_score_grid_prefix4                                 ${PDB}.resid_
multigrid_score_grid_prefix5                                 ${PDB}.resid_
multigrid_score_grid_prefix6                                 ${PDB}.resid_
multigrid_score_grid_prefix7                                 ${PDB}.resid_
multigrid_score_grid_prefix8                                 ${PDB}.resid_
multigrid_score_grid_prefix9                                 ${PDB}.resid_
multigrid_score_grid_prefix10                                ${PDB}.resid_
multigrid_score_grid_prefix11                                ${PDB}.resid_
multigrid_score_grid_prefix12                                ${PDB}.resid_
multigrid_score_grid_prefix13                                ${PDB}.resid_
multigrid_score_grid_prefix14                                ${PDB}.resid_
multigrid_score_grid_prefix15                                ${PDB}.resid_
multigrid_score_grid_prefix16                                ${PDB}.resid_
multigrid_score_grid_prefix17                                ${PDB}.resid_
multigrid_score_grid_prefix18                                ${PDB}.resid_
multigrid_score_grid_prefix19                                ${PDB}.resid_
multigrid_score_grid_prefix20                                ${PDB}.resid_
multigrid_score_grid_prefix21                                ${PDB}.resid_
multigrid_score_grid_prefix22                                ${PDB}.resid_
multigrid_score_grid_prefix23                                ${PDB}.resid_
multigrid_score_grid_prefix24                                ${PDB}.resid_
multigrid_score_grid_prefix25                                ${PDB}.resid_remaining
multigrid_score_fp_ref_mol                                   no
multigrid_score_fp_ref_text                                  yes
multigrid_score_footprint_text                               ../${PDB}.reference.txt
multigrid_score_foot_compare_type                            Euclidean
multigrid_score_normalize_foot                               no
multigrid_score_vdw_euc_scale                                1.0
multigrid_score_es_euc_scale                                 1.0
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
simplex_max_iterations                                       1000
simplex_tors_premin_iterations                               0
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                ${DOCKHOME}/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               ${DOCKHOME}/parameters/flex.defn
flex_drive_file                                              ${DOCKHOME}/parameters/flex_drive.tbl
ligand_outfile_prefix                                        output
write_footprints                                             no
write_orientations                                           no
num_scored_conformers                                        1
rank_ligands                                                 no
EOF


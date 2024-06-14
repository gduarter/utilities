#!/bin/bash

# Author: Guilherme Duarte Ramos Matos
# Date: June 2024


###### Help Function ########
helpFunction(){
    echo -e "\t\tThis script was written by Guilherme D. R. Matos based on OpenBioSim's tutorials."
    echo -e "\tIt requires Python 3.9 with the packages BioSimSpace, sire, and numpy."
    echo -e "  "
    echo -e "\tUsage: $0 -p directory_with_protein_pdb -l directory_with_ligands_sdfs -s directory_with_important_scripts"
    echo -e " "
    echo -e "\tAll directories above MUST BE in the directory where you run this script."
    exit 1
}

# Assign typed arguments to variables
while getopts ":p:l:s:" opt
do
    case $opt in
        p ) protein_dir="$OPTARG";;
        l ) ligands_dir="$OPTARG";;
        s ) scripts_dir="$OPTARG";;
        ? ) helpFunction ;;
    esac
done

# Prints helpFuntion in case the number of parameters do not match what
# the script requires
if [ -z "${protein_dir}" ] || [ -z "${ligands_dir}" ]
then
    echo "You are misusing this script"
    helpFunction
fi

python_path="/Users/guilherme/miniconda3/envs/openbiosim/bin"
root=$(pwd)
network_dir=${root}/003.network_files
md_prep_path=${root}/004.md_prep
rbfe_path=${root}/005.rbfe
scripts_path=${root}/${scripts_dir}
protein_path=${root}/${protein_dir}
ligands_path=${root}/${ligands_dir}


# Create network files
if ! [ -f ${network_dir}/network.dat ]
then
    ${python_path}/python ${scripts_path}/create_network.py -l ${ligands_dir} -n ${network_dir}
    echo "Network done!"
else
    echo " "
    echo "Network files already exist."
    echo "Check ${network_dir}"
fi

# Parameterize protein and ligands.
if ! [ -f ${md_prep_path}/md_checkpoint.out ]
then
    ${python_path}/python ${scripts_path}/prep_md_simulations.py -p ${protein_dir} -l ${ligands_dir} -o ${md_prep_path} -n ${network_dir}
    echo "MD files done!"
else
    echo " "
    echo "MD files already prepped."
    echo "Check ${md_prep_path}"
fi

# Create RBFE files.
if ! [ -f ${rbfe_path}/rbfe_checkpoint.out ]
then
    ${python_path}/python ${scripts_path}/prep_rbfe_simulations.py -n ${network_dir} -p ${md_prep_path} -r ${rbfe_path}
    echo "RBFE files done!"
else
    echo " "
    echo "RBFE files already prepped."
    echo "Check ${rbfe_path}"
fi


# Organize RBFE simulation files
types=("bound" "free")
sim_dir="lambda_0.0000" 
cd ${rbfe_path}
for dir in */
do
    for elem in "${types[@]}"
    do
        if [ -d ${rbfe_path}/${dir}/${elem}/${sim_dir} ]
        then
            cd ${rbfe_path}/${dir}/${elem}/${sim_dir}
            echo " " >> dgparameters.txt
            echo "; Free Energy Parameters" >> dgparameters.txt
            tail -n 8 gromacs.mdp >> dgparameters.txt
            sed -i 's/init-lambda-state = 0/init-lambda-state = X/' dgparameters.txt
            bash ${scripts_path}/create_gromacs_mdps_${elem}.sh
            for mdpfile in *.mdp
            do
                cat dgparameters.txt >> $mdpfile
            done
            rm gromacs.err gromacs.mdp gromacs.out gromacs.out.mdp gromacs.tpr
            mv * ${rbfe_path}/${dir}/${elem}/
            cd ..
            rm -r lambda_*
            cd ${rbfe_path}
        fi
    done
done




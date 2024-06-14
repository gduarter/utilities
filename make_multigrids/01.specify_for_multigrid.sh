#!/bin/bash

###### Help Function ########
helpFunction(){
    echo -e "\tUsage: $0 -p pdb_code -o dock6_outfile -f fp_out_txt"
    exit 1
}

# Assign typed arguments to variables
while getopts "p:o:f:" opt
do
    case $opt in
        p ) PDB="$OPTARG";;
        o ) DOCK6OUT="$OPTARG";;
        f ) FPTXT="$OPTARG";;
        ? ) helpFunction ;;
    esac
done


# Prints helpFuntion in case the number of parameters do not match what
# the script requires
if [ -z "${PDB}" ] || [ -z "${DOCK6OUT}" ] || [ -z "${FPTXT}" ]
then
    echo "You are misusing this script"
    helpFunction
fi

# Create temp file containing the numbers corresponding to the best scoring
# residues
grep -A 1 "range_union" ${DOCK6OUT} | grep -v "range_union" | grep -v "\-" | sed -e '{s/,/\n/g}' | sed -e '{s/ //g}' | sed '/^$/d' | sort -n | uniq > temp.dat

# Create list of primary residues
for i in `cat temp.dat`
do
    printf "%0*d\n" 3 $i
done > ${PDB}.primary_residues.dat

for RES in `cat temp.dat`
do
    grep "${RES} " ${FPTXT} | awk -v temp=${RES} '{if ($2 == temp) print $0;}' | awk '{print $1 "  " $3 "  " $4}' >> reference.txt
done

grep "remainder" ${FPTXT} | sed -e '{s/,/  /g}' | tr -d '\n' | awk '{print $2 "  " $3 "  " $6}' >> reference.txt 
mv reference.txt ${PDB}.reference.txt 
rm temp.dat




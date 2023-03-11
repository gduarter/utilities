import os
import sys
import pandas as pd
from copy import deepcopy


#####################
# Utility Functions #
#####################

def discover_beginning(end):
    # If you know where a molecule ends, you know when the next begin
    begin = [0] # Assume the first molecule begins at idx=0
    for idx in end:
        begin.append(idx + 1)
    # After loop is done, remove last element of list
    del begin[-1]
    return begin


def discover_ends(filename, end_pattern):
    data = open(filename, 'r+')
    lines = data.readlines()

    # Loop over all lines and find pattern
    end = []
    for idx, line in enumerate(lines):
        if end_pattern in line:
            end.append(idx)
        else: continue

    # Use 'discover_beginning' function to find
    # indices where molecules begin.
    begin = discover_beginning(end)

    return begin, end


def create_dataframe_from_mol2(filename, descriptor_list):
    data = open(filename, 'r+')
    lines = data.readlines()

    # Create empty dataframe
    df = pd.DataFrame()

    # Loop over descriptor list and lines of the mol2 file to
    # collect the data.
    for entry in descriptor_list:
        tmp = []
        for line in lines:
            if entry in line:
                if entry == 'RD_SMILES':
                    broken_line = line.split(':')
                    # Append the data of interest
                    tmp.append(broken_line[1].strip('\n'))
                else:
                    broken_line = line.split()
                    # Append the data of interest
                    tmp.append(broken_line[2])
            else: continue
        df[entry] = tmp
    col_length = len(df[descriptor_list[0]])
    #primary = filename.split("_")[3]
    #secondary = filename.split("_")[4]
    #df["primary_score"] = [primary for i in range(col_length)]
    #df["secondary_score"] = [secondary for i in range(col_length)]

    return df


##############
#    Main    #
##############

# Exit if not Python 3.X
if sys.version_info[:1] < (3,):
    sys.exit("This is a Python 3 script. Python 2.7 is deprecated and should not be used.")

# Check if input was properly given
if len(sys.argv) != 2:
    print("You should use the script in the following way:")
    print(f"{sys.argv[0]} multimol_mol2file")
    sys.exit()

# Assign user input to variable
mol2file = sys.argv[1]

# Important variables
end_pattern = '0 ROOT'
descriptors = [ "Name", "Molecular_Weight", "DOCK_Rotatable_Bonds", "Formal_Charge",
                "HBond_Acceptors", "HBond_Donors", "RD_Stereocenters", "RD_LogP", 
                "RD_SYNTHA", "RD_QED", "RD_num_of_PAINS", "RD_PAINS_names", "RD_SMILES", 
                "Continuous_Score", "desc_Continuous_es_energy", "desc_FPS_vdw_energy", 
                "desc_FPS_es_energy", "desc_FPS_num_hbond" ]
#descriptors = ["Name_DOCK","From_List","List_Rank","Name_MOE","Cluster_size",
#               "TotalScore_(FPS+DCE)","Continuous_Score","Continuous_vdw_energy",
#               "Internal_energy_repulsive","Footprint_Similarity_Score",
#               "FPS_vdw_fps","FPS_es_fps","FPS_hb_fps","FPS_vdw_fp_numres",
#               "FPS_es_fp_numres","FPS_hb_fp_numres","Num_H-bonds","DOCK_rot_bonds",
#               "Pharmacophore_Score","Property_Volume_Score","Tanimoto_Score",
#               "Hungarian_Matching_Similarity_Score","Descriptor_Score",
#               "MOE_rot_bonds","Molecular_weight","Num_chiral_centers",
#               "Lipinski_donors","Lipinski_acceptors","Lipinski_druglike",
#               "Lipinski_violations","SlogP","Formal_charge","logS",
#               "Ligand_efficiency","SMILES"]

# Create dataframe from mol2 file
data = create_dataframe_from_mol2( mol2file, descriptors )

# Save csv file
csvname = mol2file[:-5] # remove ".mol2" from the end
data.to_csv(f"{csvname}.csv", index=False)

## Separate ZINC IDs to facilitate the retrieval of specific molecules from
## the mol2 file. Write to file.
#zinc_ids = deepcopy( data["Name_DOCK"] )
#begin, end = discover_ends( mol2file, end_pattern )
#if (len(zinc_ids) == len(begin)) and (len(zinc_ids) == len(end)):
#    with open(f"positions_{csvname}.dat", "w+") as f:
#        f.write("Name\tbegin_idx\tend_idx\n")
#        for i in range(len(zinc_ids)):
#            f.write(f"{zinc_ids[i]}\t{begin[i]:6}\t{end[i]:6}\n")
#else:
#    print("There is something wrong in the mol2 file")

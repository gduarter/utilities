# Author: John Bickel
# Create Date: 09 June 2022
# Purpose: Using the RDKIT_SMILES header, find all unique molecules based on
# smiles string from an input multimol2.
# Last Edit: 09 June 2022

import argparse
import os
import sys
from operator import itemgetter

# Extracts molecules paired with their smiles strings - also adds a unique
# index to each molecule


def extract_molecule_list_pair_smiles(filename, score_type, verbosity):

    molecule_list = []
    with open(filename, "r") as f:
        found_header = False
        found_score = False
        header_start = ""
        smiles_string = ""
        score = 0
        append_val = 0
        line_list = []
        for line in f:
            # automatically takes the first ### line and will identify all
            # new molecule by that line
            if "###" in line and found_header is False:
                found_header = True
                header_start = line.split(":")[0].strip()+":"
                line_list.append(line)
                if score_type in line and found_score is False:
                    score = float(line.split(":")[1].strip())
                    found_score = True
                continue

            if header_start in line and found_header is True:
                append_val += 1  # shoddy way of tracking if at the second mol
                if found_score is False:
                    print("Unable to find the user supplied score in first "
                          "mol. Check that your score is located in molecule "
                          "header. Exiting.")
                    exit()
                if append_val >= 1:  # if we are, start appending things
                    molecule_list.append([smiles_string,
                                          line_list,
                                          append_val-2,
                                          score])
                line_list = []
                line_list.append(line)
            else:
                line_list.append(line)
            if "RD_SMILES:" in line:
                smiles_string = line.split(":")[1].strip()
            if score_type in line and found_score is False:
                score = float(line.split(":")[1].strip())
                found_score = True

        molecule_list.append([smiles_string, line_list, append_val, score])

        if (verbosity):
            print(f"Number of molecules extracted: {len(molecule_list)}.")
    return(molecule_list)


# This function takes in a list with a defined score in the header, then
# will output all of the unique molecules with the most useful score.
def unique_smiles_best_score(molecule_list, score_sign, verbosity):
    # make list of index paired with smiles string and score
    smiles_and_scores = [[i, x[0], x[3]] for i, x in enumerate(molecule_list)]
    # Sorts on SMILES string, so there's some sort of order.
    smiles_and_scores.sort(key=itemgetter(1))
    # get a unique list of smiles to do searching with
    uni_smi_list = [x[0] for x in molecule_list]
    unique_numbering = list(set(uni_smi_list))
    # Sort these unique smiles to match ordering in smiles_and_scores
    unique_numbering.sort()
    numbers_check = {}
    final_list = []
    start_val = 0

    # iterate through every unique smiles string
    for entry in unique_numbering:
        counter = 0
        numbers_check[entry] = 0
        found_mol = False
        collection_list = []  # list to keep all the molecules for each SMILES
        # Starting where we ended off on the last iteration - makes traversal
        # *much* cheaper on future iterations.
        for i in range(start_val, len(smiles_and_scores)):
            if smiles_and_scores[i][1] == entry:
                if found_mol is False:
                    found_mol = True
                start_val += 1  # new index
                numbers_check[entry] += 1
                collection_list.append(smiles_and_scores[i])
            if (smiles_and_scores[i][1] != entry and found_mol is True) or\
               (i == len(smiles_and_scores)-1):
                # Sorts all collected mols by whatever user specifies
                if score_sign == "positive":
                    collection_list.sort(key=itemgetter(2), reverse=True)
                elif score_sign == "negative":
                    collection_list.sort(key=itemgetter(2))
                # 0th entry is "best" by the given score sort parameters
                if len(collection_list) > 1:
                    for entry in collection_list:
                        print(entry)
                    print("___")
                final_list.append(collection_list[0])
                break
    # Warning in case something went wrong.
    if (len(final_list) == len(unique_numbering)) is False:
        print("WARNING: Final list does not contain all unique molecules!")

    # Collect actual mol entries based on the sorted indices.
    unique_indices = [x[0] for x in final_list]
    unique_indices.sort()
    unique_mols = []
    dupe_mols = []
    for i, mol in enumerate(molecule_list):
        if i in unique_indices:
            unique_mols.append(mol[1])
        else:
            dupe_mols.append(mol[1])

    # Some sanity checking values to ensure everything adds up.
    if (verbosity):
        print(f"Number of input molecules: {len(molecule_list)}")
        print(f"Number of unique SMILES strings: {len(unique_numbering)}")
        print(f"Number of unique indices extracted: {len(final_list)}")
        print(f"Number of unique molecules extracted: {len(unique_mols)}")
        print(f"Number of duplicate molecules determined: {len(dupe_mols)}")
        print(f"Unique + Duplicates = {len(unique_mols) + len(dupe_mols)}")
    return(unique_mols, dupe_mols)


def unique_by_smiles(molecule_list, verbosity):

    # get a unique list of smiles to double check numbers
    uni_smi_list = [x[0] for x in molecule_list]
    unique_numbering = list(set(uni_smi_list))

    list_of_seen = []
    uni_mols = []
    dupe_mols = []
    pos = 0
    # checks all molecule SMILES against the unique - makes secondary
    # list showing if they've been seen.
    for entry in molecule_list:
        if entry[0] not in list_of_seen:
            list_of_seen.append(entry[0])
            uni_mols.append(entry[1])
        else:
            dupe_mols.append(entry[1])

    if (verbosity):
        print(f"Number of unique SMILES strings: {len(unique_numbering)}")
        print(f"Number of unique molecules extracted: {len(uni_mols)}")
    return(uni_mols, dupe_mols)


def write_out_mols(molecule_list, duplicate_list, outfile, outfile_duplicates,
                   verbosity):
    with open(outfile, "w") as of:
        for entry in molecule_list:
            for line in entry:
                of.write(line)
    with open(outfile_duplicates, "w") as of:
        for entry in duplicate_list:
            for line in entry:
                of.write(line)
    if(verbosity):
        print(f"Wrote {len(molecule_list)} unique molecules to {outfile}.")
        print(f"Wrote {len(duplicate_list)} duplicate molecules to"
              f"{outfile_duplicates}.")
        print(f"This means {len(molecule_list) + len(duplicate_list)} "
              "molecules were written to file.")


if __name__ == "__main__":
    # setup all the parameters
    program_desc = "This script is used to plot a set of footprints with \
    varying types out of output."
    parser = argparse.ArgumentParser(program_desc)
    parser.add_argument("-fi", "--input_file",
                        help="Input MOL2 with RDKIT_Smiles in the header.",
                        type=str)
    parser.add_argument("-fo", "--fileout", help="Specifies the name of the "
                        "unique output file. Defaults to 'unique_out.mol2'",
                        default="unique_out.mol2", type=str)
    parser.add_argument("-fd", "--file_dupes",
                        help="Specifies the name of the duplicate output file.\
                         Defaults to 'dupe_mols.mol2'",
                        default="dupe_mols.mol2", type=str)
    parser.add_argument("-st", "--score_type",
                        help="The score as listed in the header that "
                        "will be used to choose which unique molecule to "
                        "keep. Defaults to Grid_Score.",
                        default="Grid_Score")
    parser.add_argument("-ss", "--score_sign",
                        help="The sign of the score used for comparisons "
                        "between scores. Since GRID_SCORE is default type, "
                        "this defaults to negative. (Lists sorted most"
                        "negative to most positive).",
                        choices=["positive", "negative"], default="negative")
    parser.add_argument("-v", "--verbosity",
                        help="Turns verbose output on.",
                        action="store_true")
    args = parser.parse_args()
    if (args.verbosity):
        print("Verbose output turned on. Printing to console.")
    # run the functions
    mol_list = extract_molecule_list_pair_smiles(args.input_file,
                                                 args.score_type,
                                                 args.verbosity)
    mol_list_score, dupe_list_score = unique_smiles_best_score(mol_list,
                                                               args.score_sign,
                                                               args.verbosity)
    write_out_mols(mol_list_score, dupe_list_score, args.fileout,
                   args.file_dupes, args.verbosity)

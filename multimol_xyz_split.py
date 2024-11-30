import argparse
import sys

# Classes, methods and functions
def split_list_by_indices(all_xyzlists, idx_list):
    """
    Splits a large list into sublists according to a given index list.
    
    Args:
        all_xyzlists (list): List containing all the data to be split.
        idx_list (list): Index list containing the breaking positions of the larger list.
        xyzlists (list): List of final sublists.
    """
    # Break file into sublists
    xyzlists = []
    start = 0
    for idx in idx_list:
        xyzlists.append(all_xyzlists[start:idx])
        start = idx # Update starting point
            
    # Add remaining portion of the list
    xyzlists.append(all_xyzlists[start:])
    
    return xyzlists


def split_ensemble_file(input_file, output_prefix, split_string):
    """
    Splits a multimolecule xyzfile into multiple individual molecule xyzfiles.
    It works with ORCA's GOAT-generated XYZ files.
    
    Args:
        input_file (str): Path to the input file.
        output_prefix (str): Prefix for the output files.
        split_string (str): The string that triggers a split when found in a line.
    """
    try:
        with open(input_file, 'r') as infile:
            all_xyzlists = infile.readlines()
        
        # Discover places where new files start and end
        idx_list = []
        for idx, line in enumerate(all_xyzlists):
            if split_string in line:
                idx_list.append(idx)
        idx_list = [num - 1 for num in idx_list]
        #print(idx_list) # places where lines will be broken
        idx_list = sorted(idx_list)
        if 0 in idx_list: idx_list.remove(0)
        
        # Generate sublists containing each xyz file
        xyzlists = split_list_by_indices(all_xyzlists, idx_list)
        #print(xyzlists)
        
        # Save each sublist in a different file
        for idx, xyzlist in enumerate(xyzlists):
            with open(f"{output_prefix}_{idx:03d}.xyz", "w+") as outfile:
                for line in xyzlist:
                    outfile.write(line)
                print(f"{output_prefix}_{idx:03d}.xyz created!")        
        
    except Exception as e:
        print(f"An error occurred: {e}")


# MAIN
if __name__ == '__main__':
    #define arguments read by this script
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", help="Multimolecule XYZ file",
                        type=str, nargs=1, required=True)
    parser.add_argument("-s", "--string", help="Pattern that repeats itself every line after the number of atoms in a XYZ file. If not present, this code might not be helpful. Avoid numbers as patterns.",
                        type=str, nargs=1, required=True)
    parser.add_argument("-o", "--outfile", help="Output filename",
                        type=str, nargs=1, required=True)


    #Exit if parameters were not given
    args = parser.parse_args()
    if len(args.infile) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.string) != 1:
        parser.print_help()
        sys.exit() 
    elif len(args.outfile) != 1:
        parser.print_help()
        sys.exit()


    #Get filename
    filename = args.infile[0]
    pattern = args.string[0]
    outname = args.outfile[0]

    #Generate multiple conformers
    split_ensemble_file(filename, outname, pattern)
    print("Done!")




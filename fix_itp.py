import numpy as np
import argparse
import sys, os
import shutil

# Important variables
header_size = 4

if __name__ == '__main__':

    # Fetch arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--itpfile", help="GROMACS itp file",
                        type=str, nargs=1, required=True)
    #parser.add_argument("--c_alpha", action="store_true", required=False)
    args = parser.parse_args()
    if len(args.itpfile) != 1:
        print(f"Usage: python {sys.argv[0]} -i itpfile.itp")
        parser.print_help()
        sys.exit()

    # Store header
    itpname = args.itpfile[0]
    file_lines = []
    with open(itpname, "r+") as f:
        lines = f.readlines()
        for idx in np.arange(0,header_size):
            file_lines.append(lines[idx])

    # Store matrix with itp data
    table = np.genfromtxt(itpname, usecols=(0,1,2,3,4), skip_header=header_size, dtype=int)

    # Transpose matrix to facilitate renumbering operation
    table = table.transpose()

    # Create temporary list with atom numbers
    tmp = table[0]

    # Renumber list by subtracting (tmp[0]-1) of all tmp elements
    tmp = tmp - tmp[0] + 1 # 5 if constraining alpha carbons, 1 if backbone or protein-H (gromacs files only)
    ###if -ca flag is TRUE, then tmp = tmp - tmp[0] + 5

    # Create table with updated atom numbers
    new_table = np.array([tmp, table[1], table[2], table[3], table[4]])
    table = new_table.transpose()

    # Create new list of lines
    for row in table:
        linestring = f"{row[0]: 4d}{row[1]: 5d}{row[2]: 5d}{row[3]: 5d}{row[4]: 5d}\n"
        file_lines.append(linestring)

    # Write file
    with open(f"{itpname}.tmp", "w+") as f:
        for line in file_lines:
            f.write(line)

    # Organize files
    root = os.getcwd()
    shutil.move(f"{root}/{itpname}", f"{root}/{itpname}.bckp")
    shutil.move(f"{root}/{itpname}.tmp", f"{root}/{itpname}")

    print("Fixing done")










from rdkit import Chem
import argparse
import sys

if __name__ == '__main__':
    # Define arguments read by this script
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="PDB file",
                        type=str, nargs=1, required=True)

    # Exit if parameters were not given
    args = parser.parse_args()
    if len(args.pdb) != 1:
        parser.print_help()
        sys.exit()

    # Create molecule
    filename = args.pdb[0]

    # New filename
    newName = filename[:-4]
    rdmol = Chem.MolFromPDBFile(filename, removeHs=False)
    file = Chem.SDWriter(f"{newName}.sdf")
    file.write(rdmol)
    file.close()

    print(f"File {newName}.sdf ready!")






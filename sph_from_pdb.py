import argparse
import sys

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdbfile", help="PDB file",
                        type=str, nargs=1, required=True)
    parser.add_argument("-r", "--radius", help="Sphere radii",
                        type=str, nargs=1, required=True)    
    parser.add_argument("-o", "--output", help="Output filename",
                        type=str, nargs=1, required=True)

    args = parser.parse_args()
    if len(args.pdbfile) != 1:
        print(f'''No pdbfile included.
    Usage: python {sys.argv[0]} -p pdbfile.pdb -r radius_value -o outputname''')
        parser.print_help()
        sys.exit()
    elif len(args.radius) != 1:
        print(f'''No radius value included. 
    Usage: python {sys.argv[0]} -p pdbfile.pdb -r radius_value -o outputname''')
        parser.print_help()
        sys.exit()
    elif len(args.output) != 1:
        print(f'''No output name included.
    Usage: python {sys.argv[0]} -p pdbfile.pdb -r radius_value -o outputname''')
        parser.print_help()
        sys.exit()

    pdbfile = args.pdbfile[0]
    radius = args.radius[0]
    output = args.output[0].split(".sph")[0]    

    # Read pdb file
    with open(pdbfile, "r+") as f:
        lines = f.readlines()
    
    # Define important variables
    num_atoms = len(lines)
    headline = f"cluster     1   number of spheres in cluster   {num_atoms}"
    
    # Extract coordinates from pdb file
    xcoord, ycoord, zcoord = [], [], []
    for line in lines:
        frag = line.split()
        xcoord.append(frag[6])
        ycoord.append(frag[7])
        zcoord.append(frag[8])
    
    # Create SPH file
    fill= " "
    width= 5
    with open(f"{output}.sph", "w+") as g:
        g.write(f"{headline}\n")
        for idx, tup in enumerate(zip(xcoord,ycoord,zcoord)):
            string = f'''{idx+1:{fill}{width}}{tup[0]:>10}{tup[1]:>10}{tup[2]:>10}{radius:>8}{idx+1:{fill}{width}}\n'''
            g.write(string)

print("SPH file created!")






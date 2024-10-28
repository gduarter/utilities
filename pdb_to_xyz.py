from rdkit import Chem
import argparse
import sys

class PDBtoXYZ:
    def __init__(self, pdbFile, spinMult, charge):
        self.pdbFile = pdbFile
        self.spinMult = spinMult
        self.charge = charge

    def createXYZ(self):
        rdmol = Chem.MolFromPDBFile(self.pdbFile,
                                    removeHs=False)
        numAtoms = rdmol.GetNumAtoms()
        xyzBlock = f'''{numAtoms}\n#charge={self.charge} multiplicity={self.spinMult}\n'''
        for idx, atom in enumerate(rdmol.GetAtoms()):
            pos = rdmol.GetConformer().GetAtomPosition(idx)
            xyzBlock += f'''{atom.GetSymbol():2s}  {pos.x:10.6f} {pos.y:10.6f} {pos.z:10.6f}\n'''
        return xyzBlock

    def createForPsi4(self):
        rdmol = Chem.MolFromPDBFile(self.pdbFile,
                                    removeHs=False)
        psi4Block = f'''{self.charge} {self.spinMult}\n'''
        for idx, atom in enumerate(rdmol.GetAtoms()):
            pos = rdmol.GetConformer().GetAtomPosition(idx)
            psi4Block += f'''{atom.GetSymbol():2s}  {pos.x:10.6f} {pos.y:10.6f} {pos.z:10.6f}\n'''
        return psi4Block

    def createForOrca(self):
        rdmol = Chem.MolFromPDBFile(self.pdbFile,
                                    removeHs=False)
        orcaBlock = f'''*xyz {self.charge} {self.spinMult}\n'''
        for idx, atom in enumerate(rdmol.GetAtoms()):
            pos = rdmol.GetConformer().GetAtomPosition(idx)
            orcaBlock += f'''{atom.GetSymbol():2s}  {pos.x:10.6f} {pos.y:10.6f} {pos.z:10.6f}\n'''
        orcaBlock += f'''*\n'''
        return orcaBlock


if __name__ == '__main__':
    # Define arguments read by this script
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pdb", help="PDB file",
                        type=str, nargs=1, required=True)
    parser.add_argument("-s", "--multiplicity", help="Spin multiplicity",
                        type=int, nargs=1, required=True)
    parser.add_argument("-c", "--charge", help="Molecular charge",
                        type=int, nargs=1, required=True)

    # Exit if parameters were not given
    args = parser.parse_args()
    if len(args.pdb) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.multiplicity) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.charge) != 1:
        parser.print_help()
        sys.exit()

    # Create molecule
    filename = args.pdb[0]
    multiplicity = args.multiplicity[0]
    charge = args.charge[0]

    # New filename
    newName = filename[:-4]

    # Create object and print files
    molText = PDBtoXYZ(filename, multiplicity, charge)
    f1 = molText.createXYZ()
    f2 = molText.createForPsi4()
    f3 = molText.createForOrca()
    print(f1, file = open(f'{newName}.xyz','w+'))
    print(f2, file = open(f'{newName}.psi4.xyz','w+'))
    print(f3, file = open(f'{newName}.inp','w+'))

    print("PDB turned to XYZ done!")



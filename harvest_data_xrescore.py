import os, sys
import glob
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import MACCSkeys

if sys.version_info[:1] < (3,):
    sys.exit("This is a Python 3 script. Python 2.7 is deprecated and SHOULD NOT be used")

if len(sys.argv) != 2:
    sys.exit(f"USAGE: python {sys.argv[0]} mol2file")

# Important variables
mol2file = sys.argv[1]
observables = ["RD_Stereocenters", "RD_LogP", "RD_TPSA", "RD_SYNTHA", "RD_QED", "RD_num_arom_rings",
               "RD_Spiro_atoms", "Continuous_Score", "Footprint_Similarity_Score"]

os.system("rm *_grep.txt")

# Create txt files for analysis
for elem in observables:
    try:
        os.system(f"grep '     {elem}:' {mol2file} > {elem}_grep.txt") 
    except:
        print("Check MOL2 file")

# text files containing the data
all_txt = glob.glob("*_grep.txt")

# DataFrame that will store all data
df = pd.DataFrame()
for elem in all_txt:
    tmpArr = np.genfromtxt(elem, dtype=float, usecols=2, comments="!")
    observable = elem.split("_grep")[0]
    # Fill dataframe
    df[observable] = tmpArr

# Create txt file for names
os.system(f"grep '     Name:' {mol2file} > Name_grep.txt")
names = np.genfromtxt("Name_grep.txt", dtype=str, delimiter=":", usecols=1, comments="!")
os.system("rm Name_grep.txt")

# Create txt file for SMILES
os.system(f"grep 'RD_SMILES:' {mol2file} > SMILES_grep.txt")
smi = np.genfromtxt("SMILES_grep.txt", dtype=str, delimiter=":", usecols=1, comments="!")
os.system(f"rm SMILES_grep.txt")

# Assign names to dataframe
df["Name"] = names
df["SMILES"] = smi

# Calculate MACCS keys
mols = [Chem.MolFromSmiles(elem) for elem in smi]
keys = [MACCSkeys.GenMACCSKeys(mol) for mol in mols ]
df["MACCS"] = keys

# Create name for results file
filename = mol2file.strip(".mol2")

# Save data to csv
df.to_csv(f"{filename}.csv", index=False)

# Erase grep.txt files
for elem in all_txt:
    os.system(f"rm {elem}")

print("Analysis done!")
print(" ")

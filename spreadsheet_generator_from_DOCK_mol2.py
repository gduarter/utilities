# Created by John Bickel and adapted by Guilherme Duarte

# Import modules for our exercise
import os
import sys
from PIL import Image, ImageChops
import xlsxwriter

# Import modules for dealing with chemical information
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as Chem2
from rdkit.Chem import Descriptors as Desc
from rdkit.Chem import rdmolfiles as RDFile
from rdkit.Chem import fmcs, rdFMCS, Draw, rdMolDescriptors

### Utility functions

# Function that reads multi-molecule MOL2 files. Adapted from:
# https://chem-workflows.com/articles/2019/07/18/building-a-multi-molecule-mol2-reader-for-rdkit/
# further adapted from function written by Guilherme Duarte and John Bickel, Rizzo Lab
def mol2_mol_supplier_loop(file):
    ''' This function extracts all the molecules in the multi-molecule
        MOL2 file `file` and returns a list of rdkit.Chem.rdchem.mol 
        object.
        
        Variable         I/O          dtype           default value?
        ------------------------------------------------------------
        file              I           string                  None
        mols              O           list                    N/A
        mols[i]           O           rdkit.Chem.rdchem.mol   N/A
        
    '''
    mols=[]
    smiles=[]    
    recording=False
    with open(file, 'r') as f:
        for line in f:
            
            # Determines if @<TRIPOS>MOLECULE is in line, which marks the start
            # of each molecule. Sets the recording variable to True, which is
            # the boolean to write each molecule.
            if ("@<TRIPOS>MOLECULE") in line:
                recording=True
                mol = []
            # Determines if "ROOT" is in line, which marks the end of each
            # molecule. Records the line and sets the recording variable
            # to false.
            elif ("ROOT") in line:
                mol.append(line)
                recording=False
                
                # Makes final adjustments to the data. It must look
                # like the MOL2 file of a single molecule.
                block = ",".join(mol).replace(',','')
                
                # Converts the data of a single molecule to a 
                # rdkit.Chem.rdchem.mol object.
                m=Chem.MolFromMol2Block(block,
                                        sanitize=False,
                                        cleanupSubstructures=False)
                mols.append(m)
                continue
            elif ("SMILES") in line:
                smi = line.split()[2]
                smiles.append(smi)                
                
            if recording==True:
                mol.append(line)
                
        return(mols, smiles)

##For image trimming
def trim(im):
    bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
    diff = ImageChops.difference(im, bg)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    #Bounding box given as a 4-tuple defining the left, upper, right, and lower pixel coordinates.
    #If the image is completely empty, this method returns None.
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)

###### Main

# File check
if len(sys.argv) != 2:
    print("This script requires one argument")
    print(f"Usage: {sys.argv[0]} mol2file")
    sys.exit()
filename = sys.argv[1]
outcsv_name=f"{filename[:-5]}.xlsx"


# Make directory for images, if it doesn't exist (similar to "mkdir -p")
rundir = os.getcwd()
image_dir = os.path.join(rundir, "zzz.figures")
try:
    os.makedirs(image_dir)
except: pass 

#Creates the rdkit list from an input mol2
molecule_list, smiles =mol2_mol_supplier_loop(filename)

if len(molecule_list) != len(smiles):
    print("#SMILES != #mols")
    sys.exit()

#removes hydrogens
for i in range(0,len(molecule_list)):
    try:
        molecule_list[i]=Chem.RemoveHs(molecule_list[i])
    except:
        mol = Chem.MolFromSmiles(smiles[i])
        molecule_list[i]=Chem.RemoveHs(mol)

#Computs 2d coordinates and generates a SMARTS string
for m in molecule_list: tmp=Chem2.Compute2DCoords(m)
res=rdFMCS.FindMCS(molecule_list, timeout=1).smartsString

#Legacy code for substructure matching - keeping in case it's needed
#res1=Chem.MolFromSmarts(res)
#res2=Chem2.Compute2DCoords(res1)
#subms = [x for x in molecule_list if x.HasSubstructMatch(res1)]

#Generates the image files
molecule_number=0
for m in molecule_list:
    molecule_number+=1
    
    #Legacy code for substructure matching - keeping in case it's needed
    #Chem2.GenerateDepictionMatching2DStructure(m,res1)
    
    Chem.Draw.MolToImageFile(m,'%s/%d.png' % (image_dir, molecule_number))

## Prepare figures and place them in a spreadsheet
filelist=[]
filelist_trimmed=[]
width_list=[]
height_list=[]
for i in range(1,len(molecule_list)+1):
    filelist.append("%s/%d.png" % (image_dir,i))

#print(filelist)

for i in range(0,len(filelist)):
    bg = Image.open(filelist[i]) # The image to be cropped
    #new_im = trim(bg)
    #new_im.save("%s/%s_trimmed.png" % (image_dir,str(i+1)),dpi=(96,96))
    #filelist_trimmed.append("%s/%s_trimmed.png" % (image_dir,str(i+1)))
    width_list.append(bg.width) #new_im.width)
    height_list.append(bg.height) #new_im.height)

max_width=max(width_list)
max_height=max(height_list)

#print(height_list)
#print(width_list)
#print(width_list.index(245))
divide_width=6.4
divide_height=1.4
image_row=1
image_col=5

workbook = xlsxwriter.Workbook(outcsv_name)
worksheet = workbook.add_worksheet('sheet1')

counter=1
for file in filelist: #filelist_trimmed:
    img_file = Image.open(file)
       
    # get original image parameters...
    #width, height = img_file.size
    #print(width)
    #format = img_file.format
    #mode = img_file.mode

    worksheet.set_column(image_col,image_col,(max_width+5)/divide_width)
    worksheet.set_row(image_row,(max_height+5)/divide_height)
    worksheet.insert_image(image_row, image_col, file, {'x_scale': 1, 'y_scale': 1, 
                                                        'x_offset' : 10, 'y_offset' : 10, 
                                                        'object_position':2})
    worksheet.write(image_row,image_col+1,counter)
    worksheet.write(image_row,0,counter)
    image_row+=1
    counter+=1

#worksheet.set_column(image_col+1,image_col+1,None,{'hidden': True})
#worksheet.set_default_row(hide_unused_rows=True)
#worksheet.set_column('G:XFD', None, None, {'hidden': True})
workbook.close()





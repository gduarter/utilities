import argparse, sys
from glob import glob
import os

os.environ['NUMEXPR_MAX_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

if __name__ ==  '__main__':
    import BioSimSpace as BSS
    node = BSS.Gateway.Node("A node to create input files for molecular dynamics simulation.")
    node.addInput("Ligand FF", BSS.Gateway.String( help="Force field to parameterise ligands with.", allowed=["GAFF1", "GAFF2"],default="GAFF2",),)
    node.addInput("Protein FF", BSS.Gateway.String(help="Force field to parameterise the protein with.", allowed=["FF03", "FF14SB", "FF99", "FF99SB", "FF99SBILDN"], default="FF14SB",),)
    node.addInput("Water Model", BSS.Gateway.String(help="Water model to use.", allowed=["SPC", "SPCE", "TIP3P", "TIP4P", "TIP5P"], default="TIP3P",),)
    node.addInput("Box Edges", BSS.Gateway.String(help="Size of water box around molecular system.", allowed=["20*angstrom","25*angstrom","30*angstrom","35*angstrom","45*angstrom","5*nm","7*nm","10*nm",],default="25*angstrom",),)
    node.addInput("Box Shape", BSS.Gateway.String(help="Geometric shape of water box.", allowed=["cubic", "truncatedOctahedron"], default="cubic",),)
    node.addInput("Run Time", BSS.Gateway.String(help="The sampling time per lambda window.", allowed=["10*ps","100*ps","1*ns","2*ns","3*ns","4*ns","5*ns","8*ns","10*ns","12*ns","15*ns",],default="10*ns",),)

    # add SOMD2 as an FEP engine option
    fep_engines = [e.upper() for e in BSS.FreeEnergy.engines()]
    if "SOMD2" not in fep_engines:
        fep_engines.append("SOMD2")

    node.addInput("FEP Engine", BSS.Gateway.String(help="Engine to run FEP with.", allowed=fep_engines, default="GROMACS",),)
    node.addInput("LambdaWindows", BSS.Gateway.String( help="The number of lambda windows for regular transformations.", allowed=["3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20",],default="15",),)
    node.addInput("DiffLambdaWindows", BSS.Gateway.String(help="The number of lambda windows for difficult transformations.",allowed=["4","5","6","7","8","9","10","11","12", "13","14","15","16","17", "18","19","20",],default="20",),)
    node.addInput("LOMAP Threshold", BSS.Gateway.String(help="The LOMAP score threshold to define difficult transformations.", allowed=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"],default="0.4",),)


if __name__ ==  '__main__':
    # Arguments read by this script
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--protein", help="path to protein file.",
                        type=str, nargs=1, required=True)
    parser.add_argument("-l", "--ligands", help="path to network files.",
                        type=str, nargs=1, required=True)
    parser.add_argument("-o", "--out_dir", help="path to MD directories.",
                        type=str, nargs=1, required=True)
    parser.add_argument("-n", "--network", help="path to network files.",
                        type=str, nargs=1, required=True)
    
    
    # Exit if parameters are not given
    args = parser.parse_args()
    if len(args.protein) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.ligands) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.out_dir) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.network) != 1:
        parser.print_help()
        sys.exit()
    
    # Define paths
    ligands_path = args.ligands[0]
    protein_path = args.protein[0]
    out_dir = args.out_dir[0]
    network_path = args.network[0]
    



    # Parameterize the protein
    prot = BSS.IO.readPDB(f"{protein_path}/protein.pdb", pdb4amber=False)[0]
    prot_p = BSS.Parameters.parameterise(prot, node.getInput("Protein FF")).getMolecule()
    BSS.IO.saveMolecules(f"{protein_path}/protein", prot_p, ["PRM7", "RST7"])
    print("Protein successfully parameterized.")
    
    # List all ligands
    ligand_files = sorted(glob(f"{ligands_path}/*.sdf"))
    
    # Parameterize ligands and create protein+ligand simulation files
    for idx, ligand in enumerate(ligand_files):
        lig_name = ligand.split("/")[-1][:-4]
        ### Load matching input with BSS.read.IO.
        lig = BSS.IO.readMolecules(f"{ligand}")[0]
        #################
        ### parameterise ligand by deriving requested FF from protocol.dat.
        stream = open(f"{network_path}/protocol.dat", "r")
        lines = stream.readlines()
        ff_query = lines[0].rstrip()
    
        # do tests on protocol, parameterise with correct FF.
        if not "ligand" in ff_query:
            raise NameError(
                "Please supply a ligand force field on the first line of protocol.dat. The line should look like (e.g.):\n"
                + "ligand forcefield = GAFF2"
            )
        else:
            if "GAFF1" in ff_query:
                print(f"Parameterising {lig_name} using GAFF1 force field.")
                lig_p = BSS.Parameters.gaff(lig).getMolecule()
    
            elif "GAFF2" in ff_query:
                print(f"Parameterising {lig_name} using GAFF2 force field.")
                lig_p = BSS.Parameters.gaff2(lig).getMolecule()
            else:
                raise NameError(
                    f"Force field not supported: {ff_query}. Please use either of [GAFF1, GAFF2, OpenForceField], or"
                    + " edit this script to use other force fields available in BSS.Parameters.forceFields()."
                )
        BSS.IO.saveMolecules(f"{out_dir}/ligands/{lig_name}/{lig_name}", lig_p, ["PRM7", "RST7"])
        # at this point BSS should have thrown an error if parameterisation failed, so no checks needed.
        ### make a copy of ligand, solvate original ligand object.
        lig_p_copy = lig_p.copy()
    
        # need to derive some settings from protocol.dat again.
        stream = open(f"{network_path}/protocol.dat", "r")
        lines = stream.readlines()
    
        # get the solvent force field.
        solvent_query = lines[2].rstrip().replace(" ", "").split("=")[-1]
    
        # get the box size settings.
        boxsize_query = lines[3].rstrip().replace(" ", "").split("=")[-1]
        box_axis_length = boxsize_query.split("*")[0]
        box_axis_unit = boxsize_query.split("*")[1]
        if box_axis_unit.lower() == "nm" or box_axis_unit.lower() == "nanometer":
            box_axis_unit = BSS.Units.Length.nanometer
    
        elif box_axis_unit.lower() == "a" or box_axis_unit.lower() == "angstrom":
            box_axis_unit = BSS.Units.Length.angstrom
        else:
            raise NameError(
                "Input unit not recognised. Please use any of ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p'] in "
                + "the fourth line of protocol.dat in the shape of (e.g.):\nbox edges = 10*angstrom"
            )
    
        box_min, box_max = lig_p.getAxisAlignedBoundingBox()
        box_size = [y - x for x, y in zip(box_min, box_max)]
        box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]
    
    
        # do the same for ligand +protein system.
        protein = BSS.IO.readMolecules([f"{protein_path}/protein.rst7", f"{protein_path}/protein.prm7"])[0]
        system = lig_p + protein
    
        # box size
        boxsize_query = lines[3].rstrip().replace(" ", "").split("=")[-1]
        box_axis_length = boxsize_query.split("*")[0]
        box_axis_unit = boxsize_query.split("*")[1]
        if box_axis_unit.lower() == "nm" or box_axis_unit.lower() == "nanometer":
            box_axis_unit = BSS.Units.Length.nanometer
        elif box_axis_unit.lower() == "a" or box_axis_unit.lower() == "angstrom":
            box_axis_unit = BSS.Units.Length.angstrom
        else:
            raise NameError(
                "Input unit not recognised. Please use any of ['spc', 'spce', 'tip3p', 'tip4p', 'tip5p'] in "
                + "the fourth line of protocol.dat in the shape of (e.g.):\nbox edges = 10*angstrom"
            )
    
        # get the box type settings.
        boxtype_query = lines[4].rstrip().replace(" ", "").split("=")[-1]
    
        boxtype_dict = {
            "cubic": BSS.Box.cubic,
            "truncatedoctahedron": BSS.Box.truncatedOctahedron,
            "octahedral": BSS.Box.truncatedOctahedron,
        }
    
        if boxtype_query not in boxtype_dict:
            raise NameError(
                "Input box type not recognised. Please use any of ['orthorhombic', 'octahedral', 'triclinic']"
                + "in the fifth line of protocol.dat in the shape of (e.g.):\nbox type = orthorhombic"
            )
    
    
        boxtype_func = boxtype_dict[boxtype_query]
        print(f"Solvating for {lig_name}...")
    
        box_min, box_max = lig_p.getAxisAlignedBoundingBox()
        box_size = [y - x for x, y in zip(box_min, box_max)]
        box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]
    
        box, angles = boxtype_func(max(box_sizes))
        lig_p_solvated = BSS.Solvent.solvate(
            solvent_query, molecule=lig_p, box=box, angles=angles, ion_conc=0.15
        )
        BSS.IO.saveMolecules(f"{out_dir}/solvated_ligands/{lig_name}/{lig_name}", lig_p_solvated, ["PRM7", "RST7"])
    
        box_min, box_max = system.getAxisAlignedBoundingBox()
        box_size = [y - x for x, y in zip(box_min, box_max)]
        box_sizes = [x + int(box_axis_length) * box_axis_unit for x in box_size]
    
        box, angles = boxtype_func(max(box_sizes))
        system_solvated = BSS.Solvent.solvate(
            solvent_query, molecule=system, box=box, angles=angles, ion_conc=0.15
        )
        BSS.IO.saveMolecules(f"{out_dir}/solvated_protein_ligands/{lig_name}/{lig_name}", system_solvated, ["PRM7", "RST7"])
    
        if idx == len(ligand_files) - 1:
            with open(f"{out_dir}/md_checkpoint.out", "w+") as f:
                f.write("All ligands were parameterized and protein-ligand MD files are ready")
    
    print("MD prepping done!")


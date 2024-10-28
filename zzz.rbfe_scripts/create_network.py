import os

os.environ['NUMEXPR_MAX_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

import argparse, sys, glob
import numpy as np
import pickle, csv 
if __name__ ==  '__main__':
    import BioSimSpace as BSS
    node = BSS.Gateway.Node("A node to create input files for molecular dynamics simulation.")
    node.addInput("Ligand FF", BSS.Gateway.String( help="Force field to parameterise ligands with.", allowed=["GAFF1", "GAFF2"],default="GAFF2",),)
    node.addInput("Protein FF", BSS.Gateway.String(help="Force field to parameterise the protein with.", allowed=["FF03", "FF14SB", "FF99", "FF99SB", "FF99SBILDN"], default="FF14SB",),)
    node.addInput("Water Model", BSS.Gateway.String(help="Water model to use.", allowed=["SPC", "SPCE", "TIP3P", "TIP4P", "TIP5P"], default="TIP3P",),)
    node.addInput("Box Edges", BSS.Gateway.String(help="Size of water box around molecular system.", allowed=["20*angstrom","25*angstrom","30*angstrom","35*angstrom","45*angstrom","5*nm","7*nm","10*nm",],default="45*angstrom",),)
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
    parser.add_argument("-l", "--ligands", help="path to ligands",
                        type=str, nargs=1, required=True)
    parser.add_argument("-n", "--network", help="path to network files",
                        type=str, nargs=1, required=True)
    
    # Exit if parameters are not given
    args = parser.parse_args()
    if len(args.ligands) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.network) != 1:
        parser.print_help()
        sys.exit()
    
    ligands_path = args.ligands[0]
    network_path = args.network[0]

    if not os.path.exists(f"{network_path}"):
        os.makedirs(f"{network_path}")

    # create protocol.
    protocol = [
        f"ligand forcefield = {node.getInput('Ligand FF')}",
        f"protein forcefield = {node.getInput('Protein FF')}",
        f"solvent = {node.getInput('Water Model')}",
        f"box edges = {node.getInput('Box Edges')}",
        f"box type = {node.getInput('Box Shape')}",
        f"protocol = default",
        f"sampling = {node.getInput('Run Time')}",
        f"engine = {node.getInput('FEP Engine').upper()}",]
    with open(f"{network_path}/protocol.dat", "w+") as protocol_file:
        writer = csv.writer(protocol_file)
        for prot_line in protocol:
            writer.writerow([prot_line])
    print(f"Protocol written at {network_path}/protocol.dat")

    # Create network
    ligand_files = sorted(glob.glob(f"{ligands_path}/*.sdf"))
    ligands = []
    ligand_names = []
    for filepath in ligand_files:
        # append the molecule object to a list.
        ligands.append(BSS.IO.readMolecules(filepath)[0])
        # append the molecule name to another list so that we can use the name of each molecule in our workflow.
        ligand_names.append(filepath.split("/")[-1].replace(".sdf", ""))
    
    transformations, lomap_scores = BSS.Align.generateNetwork(ligands, plot_network=False, names=ligand_names, work_dir=f"{network_path}/visualise_network",)
    
    pert_network_dict = {}
    transformations_named = [(ligand_names[transf[0]], ligand_names[transf[1]]) for transf in transformations]
    for score, transf in sorted(zip(lomap_scores, transformations_named)):
        pert_network_dict[transf] = score
    with open(f"{network_path}/pert_network_dict.pickle", "wb") as pickle_file:
        pickle.dump(pert_network_dict, pickle_file, protocol = pickle.HIGHEST_PROTOCOL)
    print("Perturbation network successfully pickled!")
    
    
    # write ligands file.
    with open(f"{network_path}/ligands.dat", "w") as ligands_file:
        writer = csv.writer(ligands_file)
        for lig in ligand_names:
            writer.writerow([lig])
    
    # write perts file. Base the lambda schedule on the file generated in the previous cell.
    np.set_printoptions(formatter={"float": "{: .4f}".format})
    
    # from protocol, derive the engine we want to use on the cluster.
    engine = node.getInput("FEP Engine").upper()
    
    with open(f"{network_path}/network.dat", "w") as network_file:
        writer = csv.writer(network_file, delimiter=" ")
    
        for pert, lomap_score in pert_network_dict.items():
            # based on the provided (at top of notebook) lambda allocations and LOMAP threshold, decide allocation.
            if lomap_score == None or lomap_score < float(node.getInput("LOMAP Threshold")):
                num_lambda = node.getInput("DiffLambdaWindows")
            else:
                num_lambda = node.getInput("LambdaWindows")
    
            # given the number of allocated lambda windows, generate an array for parsing downstream.
            lam_array_np = np.around(np.linspace(0, 1, int(num_lambda)), decimals=5)
    
            # make the array into a format readable by bash.
            lam_array = (
                str(lam_array_np)
                .replace("[ ", "")
                .replace("]", "")
                .replace("  ", ",")
                .replace("\n", "")
            )
    
            # write out both directions for this perturbation.
            writer.writerow([pert[0], pert[1], len(lam_array_np), lam_array, engine])
            writer.writerow([pert[1], pert[0], len(lam_array_np), lam_array, engine])
    
    print(f"Network successfully built!")



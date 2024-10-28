import argparse, sys
import os, csv

os.environ['NUMEXPR_MAX_THREADS'] = '1'
os.environ['NUMEXPR_NUM_THREADS'] = '1'

if __name__ ==  '__main__':
    import BioSimSpace as BSS

if __name__ ==  '__main__':
    # Arguments read by this script
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--network", help="path to network files.",
                        type=str, nargs=1, required=True)
    parser.add_argument("-p", "--moldyn", help="path to MD files.",
                        type=str, nargs=1, required=True)
    parser.add_argument("-r", "--rbfe", help="path to RBFE files.",
                        type=str, nargs=1, required=True)    
    
    # Exit if parameters are not given
    args = parser.parse_args()
    if len(args.network) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.moldyn) != 1:
        parser.print_help()
        sys.exit()
    elif len(args.rbfe) != 1:
        parser.print_help()
        sys.exit()
    
    network_path = args.network[0]
    md_path = args.moldyn[0]
    rbfe_path = args.rbfe[0]
    
    
    ### first, figure out which engine and what runtime the user has specified in protocol.
    stream = open(f"{network_path}/protocol.dat", "r")
    lines = stream.readlines()
    
    ### get the requested engine.
    engine_query = lines[7].rstrip().replace(" ", "").split("=")[-1].upper()
    if engine_query not in ["GROMACS"]:
        raise NameError(
            """Input MD engine not recognised.
            Please use 'GROMACS' on the eighth line 
            of protocol.dat in the shape of (e.g.):
            engine = GROMACS""")
    
    ### get the requested runtime.
    runtime_query = lines[6].rstrip().replace(" ", "").split("=")[-1].split("*")[0]
    try:
        runtime_query = int(runtime_query)
    except ValueError:
        raise NameError("Input runtime value not supported. Please use an integer on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")
    
    # make sure user has set ns or ps.
    runtime_unit_query = lines[6].rstrip().replace(" ", "").split("=")[-1].split("*")[1]
    if runtime_unit_query not in ["ns", "ps"]:
        raise NameError("Input runtime unit not supported. Please use 'ns' or 'ps' on the seventh line of protocol.dat in the shape of (e.g.):\nsampling = 2*ns")
    
    if runtime_unit_query == "ns":
        runtime_unit = BSS.Units.Time.nanosecond
    elif runtime_unit_query == "ps":
        runtime_unit = BSS.Units.Time.picosecond
    
    ####### Prepare simulation files
    with open(f"{network_path}/network.dat") as lambdas_file:
        reader = csv.reader(lambdas_file, delimiter=" ")
        for row in reader:
            lig_0, lig_1 = row[0], row[1]
            num_lambda = int(row[2])
            # Load equilibrated free inputs for both ligands. Complain if input not found. These systems already contain equil. waters.
            print(f"Loading ligands {lig_0} and {lig_1}.")
            ligs_path = f"{md_path}/solvated_ligands/"
            ligand_1_sys = BSS.IO.readMolecules([f"{ligs_path}{lig_0}/{lig_0}.rst7", f"{ligs_path}{lig_0}/{lig_0}.prm7",])
            ligand_2_sys = BSS.IO.readMolecules([f"{ligs_path}{lig_1}/{lig_1}.rst7", f"{ligs_path}{lig_1}/{lig_1}.prm7",])
            
            # Extract ligands.
            ligand_1 = ligand_1_sys.getMolecule(0)
            ligand_2 = ligand_2_sys.getMolecule(0)
            
            # Align ligand2 on ligand1
            print("Mapping and aligning..")
            print(ligand_1, ligand_2)
            mapping = BSS.Align.matchAtoms(ligand_1, ligand_2, complete_rings_only=True)
            inv_mapping = {v: k for k, v in mapping.items()}
            ligand_2_a = BSS.Align.rmsdAlign(ligand_2, ligand_1, inv_mapping)
            
            # Generate merged molecule.
            print("Merging..")
            merged_ligs = BSS.Align.merge(ligand_1, ligand_2_a, mapping)
            
            #### Get equilibrated waters and waterbox information for both bound and free. Get all information from lambda==0
            # Following is work around because setBox() doesn't validate correctly boxes with lengths and angles
            
            ligand_1_sys.removeMolecules(ligand_1)
            ligand_1_sys.addMolecules(merged_ligs)
            system_free = ligand_1_sys
            
            
            ################ now repeat above steps, but for the protein + ligand systems.
            # Load equilibrated bound inputs for both ligands. Complain if input not found
            print(f"Loading bound ligands {lig_0} and {lig_1}.")
            ligs_path = f"{md_path}/solvated_protein_ligands/"
            system_1 = BSS.IO.readMolecules([f"{ligs_path}{lig_0}/{lig_0}.rst7", f"{ligs_path}{lig_0}/{lig_0}.prm7",])
            system_2 = BSS.IO.readMolecules([f"{ligs_path}{lig_1}/{lig_1}.rst7", f"{ligs_path}{lig_1}/{lig_1}.prm7",])
            
            # Extract ligands and protein. Do this based on nAtoms and nResidues, as sometimes
            # the order of molecules is switched, so we can't use index alone.
            # bugfix in BSS makes the below redundant but keeping this in to be 100% sure we're getting the correct structures.
            system_ligand_1 = None
            protein = None
            n_residues = [mol.nResidues() for mol in system_1]
            n_atoms = [mol.nAtoms() for mol in system_1]
            for i, (n_resi, n_at) in enumerate(zip(n_residues[:20], n_atoms[:20])):
                if n_resi == 1 and n_at > 5:
                    system_ligand_1 = system_1.getMolecule(i)
                elif n_resi > 1:
                    protein = system_1.getMolecule(i)
                else:
                    pass
            
            # loop over molecules in system to extract the ligand
            system_ligand_2 = None
            
            n_residues = [mol.nResidues() for mol in system_2]
            n_atoms = [mol.nAtoms() for mol in system_2]
            for i, (n_resi, n_at) in enumerate(zip(n_residues, n_atoms)):
                # grab the system's ligand and the protein. ignore the waters.
                if n_resi == 1 and n_at > 5:
                    system_ligand_2 = system_2.getMolecule(i)
                else:
                    pass
            
            # extract ions.
            # ions_bound = system_2.search("not mols with atomidx 2")
            
            if system_ligand_1 and system_ligand_2 and protein:
                print("Using molecules ligand_1, ligand_2, protein:")
                print(system_ligand_1, system_ligand_2, protein)
            else:
                raise _Exceptions.AlignmentError("Could not extract ligands or protein from input systems. Check that your ligands/proteins are properly prepared by BSSligprep.sh!")
            
            # Align ligand2 on ligand1
            print("Mapping..")
            mapping = BSS.Align.matchAtoms(system_ligand_1, system_ligand_2, complete_rings_only=True)
            inv_mapping = {v: k for k, v in mapping.items()}
            
            print("Aligning..")
            system_ligand_2_a = BSS.Align.rmsdAlign(system_ligand_2, system_ligand_1, inv_mapping)
            
            # Generate merged molecule.
            print("Merging..")
            system_merged_ligs = BSS.Align.merge(system_ligand_1, system_ligand_2_a, mapping)
            
            system_1.removeMolecules(system_ligand_1)
            system_1.addMolecules(system_merged_ligs)
            system_bound = system_1
    
            # define the free energy protocol with all this information. User could customise settings further here, see docs.
            freenrg_protocol = BSS.Protocol.FreeEnergy(
                num_lam=num_lambda, runtime=runtime_query * runtime_unit
            )
    
            # for GROMACS, also need per lambda min + eq
            if engine_query == "GROMACS":
                min_protocol = BSS.Protocol.FreeEnergyMinimisation(num_lam=num_lambda)
            #     #nvt_protocol = BSS.Protocol.FreeEnergyEquilibration(
            #     #    num_lam=num_lambda, pressure=None, temperature= 298.15 * BSS.Units.Pressure.kelvin
            #     #)
            #     #npt_protocol = BSS.Protocol.FreeEnergyEquilibration(
            #     #    num_lam=num_lambda, pressure=1 * BSS.Units.Pressure.atm, temperature= 298.15 * BSS.Units.Pressure.kelvin
                #)
            
            ############# Set up the directory environment.
            workdir = f"{rbfe_path}/{lig_0}~{lig_1}/"
            print(f"Setting up {engine_query} directory environment in {workdir}.")
            
            print("Bound..")
            workdir = f"{rbfe_path}/{lig_0}~{lig_1}"
            BSS.FreeEnergy.Relative(
                system_bound,
                min_protocol,
                engine=f"{engine_query}",
                work_dir=workdir + "/bound/",
            )
        
            print("Free..")
            workdir = f"{rbfe_path}/{lig_0}~{lig_1}"
            BSS.FreeEnergy.Relative(
                system_free,
                min_protocol,
                engine=f"{engine_query}",
                work_dir=workdir + "/free/",
            )

    with open(f"{rbfe_path}/rbfe_checkpoint.out", "w+") as f:
        f.write(f"RBFE files are ready and can be found at {rbfe_path}") 

from pdb_reader import Protein
import numpy as np
import os
import argparse
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

def get_residue_name(prot, residue_id):
    for model in prot.structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == residue_id:
                    return residue.resname
    return None

def get_residues(prot):
    coordinates = []
    atom_names = []                
    chains = []
    residue_ids = []
    for model in prot.structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chains.append(chain.id)
                    atom_names.append(atom.get_name())
                    coordinates.append(atom.coord)
                    residue_ids.append(residue.id[1])

    def hydrogen_list():
        for i in range(len(atom_names)):
            first_letter = atom_names[i][0]
            if first_letter == "H":
                atom_names[i] = "hydrogen"
            elif first_letter.isnumeric():
                if atom_names[i][1] == "H":
                    atom_names[i] = "hydrogen"
                else:
                    atom_names[i] == "not_hydrogen"
            else:
                atom_names[i] = "not_hydrogen"

    hydrogen_list()

    listA = []
    listB = []
    residue_ids_A = []
    residue_ids_B = []
    lst = prot.get_chain_ids()
    for i in range(len(coordinates)):
        if ((chains[i] == lst[0]) and (atom_names[i] == "not_hydrogen")):
            listA.append(coordinates[i])
            residue_ids_A.append(residue_ids[i])
        elif ((chains[i] == lst[1]) and (atom_names[i] == "not_hydrogen")):
            listB.append(coordinates[i])
            residue_ids_B.append(residue_ids[i])      

    dist_thresh = 5
    res_in_contact = []
    for i,posA in enumerate(listA):
        for j,posB in enumerate(listB):
            if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                cont = [residue_ids_A[i], residue_ids_B[j]]
                if cont not in res_in_contact:
                    res_in_contact.append(cont)

    return res_in_contact


def process_pdb_folder(full_folder_path, pdb_id, hdf5_file):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    pdb_group = hdf5_file.create_group(pdb_id)
    
    for path in paths:
        print(f"Path is {path}")
        for filename in os.listdir(path):
            print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                print(f"Processing {filename}")
                protein = Protein(pdb_path)
                res_list = get_residues(protein)
                resA=[]
                resB=[]
                for res in res_list:
                    resA.append(res[0])
                    resB.append(res[1])
                for i in range(len(resA)):
                    resnameA = get_residue_name(protein, resA[i])
                    resnameB = get_residue_name(protein, resB[i])
                    results.append((resnameA, resnameB))

                model_group = pdb_group.create_group(filename)
                model_group.create_dataset('res_contacts', data=np.array(results, dtype=h5py.string_dtype()))

            else:
                print(f"File did not pass requirements.")


        

def main(folder_path):
    main_folder_name = os.path.basename(os.path.normpath(folder_path))
    hdf5_filename = f"{main_folder_name}_all_contacts.hdf5"
    print(f"CREATED FILE: {hdf5_filename}")
    hdf5_filepath = os.path.join('/vast/palmer/scratch/ohern/sr2562/hydro_results/contacts', hdf5_filename)

    with h5py.File(hdf5_filepath, 'a') as hdf5_file:
        for folder in os.listdir(folder_path):
            full_folder_path = os.path.join(folder_path, folder)
            if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
                pdb_id = full_folder_path.split('_')[-1]
                print(f"PDB id is {pdb_id}")
                process_pdb_folder(full_folder_path, pdb_id, hdf5_file)
                print("DONE----------------------------------------------------------------------")



if __name__ == "__main__":
    main(args.pdb_folder)

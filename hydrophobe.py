from pdb_reader import Protein
import numpy as np
import os
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    atom_names = []
    chains = []
    residue_ids = []
    coordinates = []

    for model in prot.structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chains.append(chain.id)
                    atom_names.append(atom.get_name())
                    coordinates.append(atom.coord)
                    residue_ids.append(residue.id[1])

    atom_names = np.array(atom_names)
    coordinates = np.array(coordinates)
    chains = np.array(chains)
    residue_ids = np.array(residue_ids)

    mask = np.char.startswith(atom_names, 'H') | np.char.startswith(atom_names, '1H') | np.char.startswith(atom_names, '2H') | np.char.startswith(atom_names, '3H')
    atom_names[mask] = 'hydrogen'
    atom_names[~mask] = 'not_hydrogen'

    listA, listB, residue_ids_A, residue_ids_B = [], [], [], []
    lst = prot.get_chain_ids()
    maskA = (chains == lst[0]) & (atom_names == 'not_hydrogen')
    maskB = (chains == lst[1]) & (atom_names == 'not_hydrogen')

    listA = coordinates[maskA]
    residue_ids_A = residue_ids[maskA]
    listB = coordinates[maskB]
    residue_ids_B = residue_ids[maskB]

    dist_thresh = 5
    res_in_contact = []

    for i, posA in enumerate(listA):
        distances = np.linalg.norm(listB - posA, axis=1)
        indices = np.where(distances <= dist_thresh)[0]
        for idx in indices:
            res_in_contact.append([residue_ids_A[i], residue_ids_B[idx]])

    return res_in_contact


def get_hydro_hits(file):
    hydrophobicity_dict = {
        "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
        "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
        "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
        "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
    }

    hits = 0

    protein = Protein(file)
    res_list = get_residues(protein)
    resA = [res[0] for res in res_list]
    resB = [res[1] for res in res_list]

    if len(resA) == len(resB):
        for i in range(len(resA)):
            resnameA = get_residue_name(protein, resA[i])
            resnameB = get_residue_name(protein, resB[i])
            if abs(hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]) <= 0.3:
                hits += 1

    return hits


def process_pdb_file(pdb_path):
    filename = os.path.basename(pdb_path)
    result = (filename[:-4], get_hydro_hits(pdb_path))
    print(result)
    return result


def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    pdb_files = []

    for path in [relaxed_folder_path, random_folder_path]:
        for filename in os.listdir(path):
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_files.append(os.path.join(path, filename))

    with ProcessPoolExecutor() as executor:
        future_to_pdb = {executor.submit(process_pdb_file, pdb_file): pdb_file for pdb_file in pdb_files}
        for future in as_completed(future_to_pdb):
            try:
                result = future.result()
                results.append(result)
            except Exception as exc:
                print(f'Generated an exception: {exc}')

    output_csv = f'{pdb_id}_hydrophobicity_contacts.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'hydrophobicity_contacts'])
        for result in results:
            writer.writerow(result)


def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-4:]
            process_pdb_folder(full_folder_path, pdb_id)

if __name__ == "__main__":
    main(args.pdb_folder)
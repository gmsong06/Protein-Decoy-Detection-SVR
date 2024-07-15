import os
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.spatial import KDTree
import multiprocessing as mp
import argparse
import csv

# Define the hydrophobicity dictionary
hydrophobicity_dict = {
    "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
    "GLN": 0.72278272, "PRO": 0.65123555, "HIS": 0.48907553, "SER": 0.52365422, "THR": 0.47798833,
    "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
    "TRP": 0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
}

# Define a function to extract residues and their coordinates from a structure
def get_residues_and_coords(structure):
    residues = []
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.id[1]
                res_name = residue.resname
                res_coords = [atom.coord for atom in residue if atom.element != 'H']
                if res_coords:
                    residues.append((chain.id, res_id, res_name))
                    coords.append(res_coords)
    return residues, coords

# Define a function to get residue pairs in contact
def get_contact_residues(coords_A, coords_B, residues_A, residues_B, dist_thresh=5.0):
    kdtree_A = KDTree(np.vstack(coords_A))
    kdtree_B = KDTree(np.vstack(coords_B))
    contact_residues_A = set()
    contact_residues_B = set()
    for i, res_coords_A in enumerate(coords_A):
        for j, res_coords_B in enumerate(coords_B):
            if kdtree_A.query_ball_tree(kdtree_B, dist_thresh):
                contact_residues_A.add(residues_A[i][1])
                contact_residues_B.add(residues_B[j][1])
    return list(contact_residues_A), list(contact_residues_B)

# Define a function to calculate hydrophobic hits
def get_hydro_hits(prot):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(prot, prot)
    residues, coords = get_residues_and_coords(structure)

    chain_A_residues = [res for res in residues if res[0] == 'A']
    chain_B_residues = [res for res in residues if res[0] == 'B']
    chain_A_coords = [coords[i] for i, res in enumerate(residues) if res[0] == 'A']
    chain_B_coords = [coords[i] for i, res in enumerate(residues) if res[0] == 'B']

    contacts_A, contacts_B = get_contact_residues(chain_A_coords, chain_B_coords, chain_A_residues, chain_B_residues)

    hits = 0
    for res_id_A, res_id_B in zip(contacts_A, contacts_B):
        res_name_A = next(res[2] for res in residues if res[1] == res_id_A and res[0] == 'A')
        res_name_B = next(res[2] for res in residues if res[1] == res_id_B and res[0] == 'B')
        if abs(hydrophobicity_dict[res_name_A] - hydrophobicity_dict[res_name_B]) <= 0.3:
            hits += 1

    return hits

# Define a function to process each PDB folder
def process_pdb_folder(args):
    print (args)
    folder_path, pdb_id = args
    results = []
    relaxed_folder_path = os.path.join(folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    for path in [relaxed_folder_path, random_folder_path]:
        for filename in os.listdir(path):
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                hits = get_hydro_hits(pdb_path)
                results.append((filename[:-4], hits))

    output_csv = f'{pdb_id}_hydrophobicity_contacts.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'hydrophobicity_contacts'])
        for result in results:
            writer.writerow(result)

# Main function to process all folders in parallel
def main(folder_path):
    tasks = []
    for folder in os.listdir(folder_path):
        if folder.startswith("sampled_") and os.path.isdir(os.path.join(folder_path, folder)):
            pdb_id = folder.split('_')[-1]
            #print(f"PATH IS {os.path.join(folder_path, folder)}")
            tasks.append((os.path.join(folder_path, folder), pdb_id))

    with mp.Pool(mp.cpu_count()) as pool:
        #print(f"TASKS IS {tasks}")
        pool.map(process_pdb_folder, tasks)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
    args = parser.parse_args()
    main(args.pdb_folder)

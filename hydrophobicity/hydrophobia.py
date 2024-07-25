from pdb_reader import Protein
import math
import numpy as np
import pandas as pd
import Bio
from Bio.PDB import PDBParser,vectors
from pathlib import Path
import os
import argparse
import csv
from collections import deque

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
    resA = []
    resB = []
    print(len(res_list))
    for res in res_list:
        resA.append(res[0])
        resB.append(res[1])

    print(f"Contacts in Chain A: {len(resA)} + Contacts in Chain B: {len(resB)}")
    if len(resA) == len(resB):
        print("Lengths are equal. Proceeding to compare hydrophobicities.")
    print("-----------------------------------------------------")

    for i in range(len(resA)):
        resnameA = get_residue_name(protein, resA[i])
        resnameB = get_residue_name(protein, resB[i])
        print(f"Contact between {resnameA} and {resnameB}")
        print(f"Hydrophobicities = {hydrophobicity_dict[resnameA]} and {hydrophobicity_dict[resnameB]}")
        if abs(hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]) <= 0.4:
            hits += 1

    return hits


def bfs(adj_list, start_node, hydro_reaction, visited, distances):
    q = deque([start_node])
    visited[start_node] = True
    distances[start_node] = 0

    while q:
        current_node = q.popleft()

        for neighbor in adj_list[current_node]:
            if not visited[neighbor]:
                visited[neighbor] = True
                distances[neighbor] = distances[current_node] + 1
                q.append(neighbor)

def get_max_dist(adj_list, hydro_reaction):
    max_distance = 0
    nodes_with_reaction = [node for node in hydro_reaction if hydro_reaction[node]]
    
    for start_node in nodes_with_reaction:
        # Initialize visited and distance dictionaries
        visited = {node: False for node in adj_list}
        distances = {node: float('inf') for node in adj_list}
        
        bfs(adj_list, start_node, hydro_reaction, visited, distances)
        
        # Calculate maximum distance to nodes with hydro_reaction as True
        for end_node in nodes_with_reaction:
            if distances[end_node] != float('inf'):
                max_distance = max(max_distance, distances[end_node])
    
    return max_distance



def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    for path in paths:
        print(f"Path is {path}")
        for filename in os.listdir(path):
            print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                print(f"Processing {filename}")
                results.append((filename[:-4], (get_hydro_hits(pdb_path))))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'/vast/palmer/scratch/ohern/sr2562/hydro_results/D/{pdb_id}_hydrophobicity_contacts.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'hydrophobicity_contacts'])
        for result in results:
            writer.writerow(result)

def main(folder_path):
    # for folder in os.listdir(folder_path):
    #     full_folder_path = os.path.join(folder_path, folder)
    #     if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
    #         pdb_id = full_folder_path[-4:]
    #         print(f"PDB id is {pdb_id}")
    #         process_pdb_folder(full_folder_path, pdb_id)
    #         print("DONE----------------------------------------------------------------------")

    graph = {
        0: [1, 4],
        1: [0, 2, 4],
        2: [1, 3, 4],
        3: [2, 5],
        4: [0, 1, 2],
        5:  [3]
    }

    hydro_reaction = {
        0: True, 
        1: True, 
        2: False,
        3: True,
        4: True,
        5: False
    }
    
    print(get_max_dist(graph, hydro_reaction))

if __name__ == "__main__":
    main(args.pdb_folder)

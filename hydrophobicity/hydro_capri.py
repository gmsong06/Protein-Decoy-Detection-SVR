from ..pdb_reader import Protein
import math
import numpy as np
import pandas as pd
import Bio
from Bio.PDB import PDBParser,vectors
from pathlib import Path
import os
import argparse
import csv

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
        if abs(hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]) <= 0.2:
            hits += 1

    return hits


def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"sampled_{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/random_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    for path in paths:
        print(f"Path is {path}")
        for filename in os.listdir(path):
            print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                print(f"Processing {filename}")
                results.append((filename[:-3], (get_hydro_hits(pdb_path))))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'/vast/palmer/scratch/ohern/sr2562/{pdb_id}_hydrophobicity_capri.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'hydrophobicity_contacts'])
        for result in results:
            writer.writerow(result)

def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-3:]
            print(f"PDB id is {pdb_id}")
            process_pdb_folder(full_folder_path, pdb_id)
            print("DONE----------------------------------------------------------------------")

if __name__ == "__main__":
    main(args.pdb_folder)
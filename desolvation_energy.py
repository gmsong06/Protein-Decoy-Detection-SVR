

"""
Protein1 = Protein("targets/1acb_complex_H.pdb")
print(Protein1.contact_dict())

def contact_dict():
    coordinates = []
    atom_ids = []
    atom_names = []
    chains = []
    residue_ids = []

    for model in self.structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chains.append(chain.id)
                    atom_ids.append(atom.get_serial_number())  # Collect atom IDs
                    atom_names.append(atom.get_name())
                    coordinates.append(atom.coord)
                    residue_ids.append(residue.id[1])

    def classify_atoms():
        for i in range(len(atom_names)):
            first_letter = atom_names[i][0]
            if first_letter == "H":
                atom_names[i] = "hydrogen"
            elif first_letter.isnumeric() and atom_names[i][1] == "H":
                atom_names[i] = "hydrogen"
            else:
                atom_names[i] = "not_hydrogen"

    classify_atoms()

    chain_dict = defaultdict(list)
    atom_dict = defaultdict(list)
    residue_dict = defaultdict(list)

    for i in range(len(coordinates)):
        if atom_names[i] == "not_hydrogen":
            chain_dict[chains[i]].append(coordinates[i])
            atom_dict[chains[i]].append(atom_ids[i])
            residue_dict[chains[i]].append(residue_ids[i])

    chain_ids = list(chain_dict.keys())
    if len(chain_ids) < 2:
        return 0

    listA = np.array(chain_dict[chain_ids[0]])
    listB = np.array(chain_dict[chain_ids[1]])

    atomsA = atom_dict[chain_ids[0]]
    atomsB = atom_dict[chain_ids[1]]

    residuesA = residue_dict[chain_ids[0]]
    residuesB = residue_dict[chain_ids[1]]

    # Using KDTree for efficient distance calculations
    treeA = KDTree(listA)
    dist_thresh = 5
    atom_residue_dict = {}

    for i, posB in enumerate(listB):
        indices = treeA.query_ball_point(posB, dist_thresh)
        for idx in indices:
            atom_residue_dict[atomsA[idx]] = residuesA[idx]
            atom_residue_dict[atomsB[i]] = residuesB[i]

    return atom_residue_dict





if residue_id = 

"""

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


def get_electro_hits(file):

    electrostatic_dict = {
            "ARG": -157.49, "ASP": -119.47, "GLU": -112.74, "LYS": -132.27, "ASN": -73.62,
            "GLN": -82.02, "PRO": -75.18, "HIS":  -96.46, "SER": -77.15, "THR": -74.43,
            "GLY": -82.95, "TYR": -73.4, "ALA": -77.40, "CYS": -74.42, "MET": -72.85,
            "TRP":  -74.11, "VAL": -71.3, "PHE": -72.71, "LEU": -81.13, "ILE": -69.51
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
        print("Lengths are equal. Proceeding to compare electrostatics.")
    print("-----------------------------------------------------")

    for i in range(len(resA)):
        resnameA = get_residue_name(protein, resA[i])
        resnameB = get_residue_name(protein, resB[i])
        print(f"Contact between {resnameA} and {resnameB}")
        print(f"Electrostatics = {electrostatic_dict[resnameA]} and {electrostatic_dict[resnameB]}")
        if abs(electrostatic_dict[resnameA] - electrostatic_dict[resnameB]) <= 45:
            hits += 1

    return hits


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
                results.append((filename[:-4], (get_electro_hits(pdb_path))))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'{pdb_id}electrostatics_contacts.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'electrostatics_contacts'])
        for result in results:
            writer.writerow(result)

def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-4:]
            print(f"PDB id is {pdb_id}")
            process_pdb_folder(full_folder_path, pdb_id)
            print("DONE----------------------------------------------------------------------")

if __name__ == "__main__":
    main(args.pdb_folder)
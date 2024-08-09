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


def get_charge_hits(file):

    charges_dict = {
            "ARG": 1, "ASP": -1, "GLU": -1, "LYS": 1, "ASN": 0,
            "GLN": 0, "PRO": 0, "HIS":  .1, "SER": 0, "THR": 0,
            "GLY": 0, "TYR": 0, "ALA": 0, "CYS": 0, "MET": 0,
            "TRP":  0, "VAL": 0, "PHE": 0, "LEU": 0, "ILE": 0
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

    if len(resA) == len(resB):
        print("Comparing Residue A charges to Residue B charges")
    print("-----------------------------------------------------")

    for i in range(len(resA)):
        resnameA = get_residue_name(protein, resA[i])
        resnameB = get_residue_name(protein, resB[i])
        print(f"Contact between {resnameA} and {resnameB}")
        print(f"Charges = {charges_dict[resnameA]} and {charges_dict[resnameB]}")
        if charges_dict[resnameA] == 0:
            hits += 0
        elif charges_dict[resnameB] == 0:
            hits += 0
        elif charges_dict[resnameA] > .1 and charges_dict[resnameB] > .1:
            hits += -1
        elif charges_dict[resnameA] < .1 and charges_dict[resnameB] < .1:
            hits += -1
        elif charges_dict[resnameA] > .1 and charges_dict[resnameB] < .1:
            hits += 1
        elif charges_dict[resnameA] < .1 and charges_dict[resnameB] > .1:
            hits += 1
        elif charges_dict[resnameA] == .1 and charges_dict[resnameB] < 0:
            hits += .1
        elif charges_dict[resnameA] == .1 and charges_dict[resnameB] > 0:
            hits += -.1
        elif charges_dict[resnameA] > 0 and charges_dict[resnameB] == .1:
            hits += -.1
        elif charges_dict[resnameA] < 0 and charges_dict[resnameB] == .1:
            hits += .1

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
                results.append((filename[:-4], (get_charge_hits(pdb_path))))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'{pdb_id}charges_score.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'charges_score'])
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
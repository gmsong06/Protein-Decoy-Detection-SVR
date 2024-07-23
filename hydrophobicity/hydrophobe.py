import numpy as np
import os
import argparse
import csv
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

def process_hdf5_file(pdb_id, model_id, hdf5_filepath):
    results = []

    hydrophobicity_dict = {
            "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
            "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
            "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
            "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
        }

    with h5py.File(hdf5_filepath, 'r') as hdf5_file:
        pdb_group = hdf5_file[pdb_id]
        model_group = pdb_group[model_id]
        res_contacts = model_group['res_contacts'][:]
        print(f"Processing Model: {model_id}")

        for contact in res_contacts:
            resnameA, resnameB = contact
            # print(f"Contact between {resnameA} and {resnameB}")
            # print(f"Hydrophobicities = {hydrophobicity_dict[resnameA]} and {hydrophobicity_dict[resnameB]}")
            fx = abs((hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]))
            results.append(fx)

    return np.mean(results) if results else None


def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    hdf5_filepath = f'/vast/palmer/scratch/ohern/sr2562/hydro_results/D/{pdb_id}_hydrophobicity_contacts.hdf5'

    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    for path in paths:
        print(f"Path is {path}")
        for filename in os.listdir(path):
            print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                model_id = filename[:-4]
                print(f"Processing {filename}")
                avg_fx = process_hdf5_file(pdb_id, model_id, hdf5_filepath)
                results.append((model_id, avg_fx))
            else:
                print(f"File did not pass requirements.")

    output_csv = f'/vast/palmer/scratch/ohern/sr2562/hydro_results/{fnc}/{pdb_id}_hydrophobicity_fnc.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'hydrophobicity_fnc'])
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
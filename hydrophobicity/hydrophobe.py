import numpy as np
import os
import argparse
import csv
import h5py
import math

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

def process_hdf5_file(hdf5_filepath):
    hydrophobicity_dict = {
            "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
            "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
            "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
            "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
        }
    
    with h5py.File(hdf5_filepath, 'r') as f:
        pdb_ids = list(f.keys())
        for pdb_id in pdb_ids:
            results = []
            models = list(f[pdb_id].keys())
            for model in models:
                res_contacts = f[f"{pdb_id}/{model}/res_contacts"][:]
                lst = []
                decoded_contacts = [[resA.decode('utf-8'), resB.decode('utf-8')] for resA, resB in res_contacts]

                for resnameA, resnameB in decoded_contacts:
                    # print(f"Contact between {resnameA} and {resnameB}")
                    # print(f"Hydrophobicities = {hydrophobicity_dict[resnameA]} and {hydrophobicity_dict[resnameB]}")
                    fx = abs((hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]))
                    lst.append(fx)

                fx_mean = np.mean(lst) if lst else None
                my_sum = 0
                for item in lst:
                    my_sum += (item - fx_mean)**2
            
                results.append([model[:-4], fx_mean, len(lst), math.sqrt((my_sum/(len(lst)-1)))])

        output_csv = f'/vast/palmer/scratch/ohern/sr2562/hydro_results/SD/{pdb_id}_hydrophobicity_fnc.csv'
        with open(output_csv, mode='w', newline='') as file:

            writer = csv.writer(file)
            writer.writerow(['pdb_file', 'avg_dif', 'N', 'standard_deviation'])
            for result in results:
                writer.writerow(result)



def main(folder_path):
    for file in os.listdir(folder_path):
        hdf5_filepath = os.path.join(folder_path, file)
        print(f"Processing {file}")
        process_hdf5_file(hdf5_filepath)
        print("DONE----------------------------------------------------------------------")

if __name__ == "__main__":
    main(args.pdb_folder)
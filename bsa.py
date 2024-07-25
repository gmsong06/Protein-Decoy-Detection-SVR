from pdb_reader import Protein
import math
import numpy as np
import pandas as pd
import Bio
from Bio.PDB import PDBParser, vectors
from pathlib import Path
import os
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

def get_bsa(pdb_path: str):
    protein = Protein(pdb_path)

    return protein.get_interface_sa()


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
                results.append((filename[:-4], get_bsa(pdb_path)))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/all_csvs/bsa/{pdb_id}_bsa.csv'
    with open(output_csv, mode='w', newline='') as file:

        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'bsa'])
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

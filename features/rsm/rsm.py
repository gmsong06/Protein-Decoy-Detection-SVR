import argparse
import freesasa
import os
import csv

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folders", type=str, help="Path to the folder containing folders of PDB files.")
args = parser.parse_args()

def calculate_Fa(SASAUnbound, SASABound):
    return ((sum(SASAUnbound)) - SASABound) / (sum(SASAUnbound))

def findUnboundSASA(pdb_path: str):
    structures = freesasa.structureArray(pdb_path, {'separate-chains': True})
    areas = []
    for structure in structures:
        result = freesasa.calc(structure)
        areas.append(result.totalArea())
    return areas

def findBoundSASA(pdb_path: str):
    structure = freesasa.Structure(pdb_path)
    result = freesasa.calc(structure)
    return result.totalArea()

def process_file(pdb_path):
    bound = findBoundSASA(pdb_path)
    unbound = findUnboundSASA(pdb_path)
    return calculate_Fa(unbound, bound)

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
                pdb_path = os.path.join(path, filename)  # Use path to get the correct file path
                fa_value = process_file(pdb_path)
                print(f"Processing {filename}")
                results.append((filename[:-4], fa_value))
            else:
                print(f"File did not pass requirements.")
    output_csv = f'{pdb_id}_interface_rsm.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'interface_rsm'])
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
    main(args.pdb_folders)

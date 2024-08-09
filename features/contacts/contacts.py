from Bio.PDB import PDBParser
import numpy as np
import argparse
import csv
import os
from multiprocessing import Pool, cpu_count
from scipy.spatial import KDTree

parser_1 = argparse.ArgumentParser()
parser_1.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser_1.parse_args()

def distance_checker(pdb_file):
    coordinates = []
    atom_names = []
    chains = []
    residue_ids = []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)

    for model in structure:
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

    chain_dict = {}
    residue_dict = {}

    for i in range(len(coordinates)):
        if chains[i] not in chain_dict:
            chain_dict[chains[i]] = []
        if chains[i] not in residue_dict:
            residue_dict[chains[i]] = []

        if atom_names[i] == "not_hydrogen":
            chain_dict[chains[i]].append(coordinates[i])
            residue_dict[chains[i]].append(residue_ids[i])

    chain_ids = list(chain_dict.keys())
    if len(chain_ids) < 2:
        return 0

    listA = np.array(chain_dict[chain_ids[0]])
    listB = np.array(chain_dict[chain_ids[1]])

    residueA = residue_dict[chain_ids[0]]
    residueB = residue_dict[chain_ids[1]]

    # Using KDTree for efficient distance calculations
    treeA = KDTree(listA)
    dist_thresh = 5
    res_in_contact = []

    for i, posB in enumerate(listB):
        indices = treeA.query_ball_point(posB, dist_thresh)
        for idx in indices:
            cont = [residueA[idx], residueB[i]]
            res_in_contact.append(cont)

    contacts = len(res_in_contact)
    return contacts

def process_pdb_file(args):
    path, filename = args
    if filename.endswith('.pdb') and ("NoH" not in filename):
        pdb_path = os.path.join(path, filename)
        return filename[:-4], distance_checker(pdb_path)
    return None

def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    pool = Pool(cpu_count())
    tasks = [(path, filename) for path in paths for filename in os.listdir(path)]
    results = pool.map(process_pdb_file, tasks)
    pool.close()
    pool.join()

    output_csv = f'{pdb_id}_all_contacts.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'all_contacts'])
        for result in results:
            if result:
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

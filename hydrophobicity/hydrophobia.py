from pdb_reader import Protein
import numpy as np
import pandas as pd
import os
import argparse
import csv

# parser = argparse.ArgumentParser()
# parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
# args = parser.parse_args()


def get_residue_name(protein, residue_id):
    for model in protein.structure:
        for chain in model:
            for residue in chain:
                if residue.id[1] == residue_id:
                    return residue.resname
    return None

def external_contacts(protein):
    coordinates = []
    atom_names = []
    chains = []
    residue_ids = []

    for model in protein.structure:
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
                    atom_names[i] = "not_hydrogen"
            else:
                atom_names[i] = "not_hydrogen"

    hydrogen_list()

    dist_thresh = 5  # Distance threshold for contacts
    
    listA = []
    listB = []
    residue_ids_A = []
    residue_ids_B = []
    lst = protein.get_chain_ids()
    for i in range(len(coordinates)):
        if ((chains[i] == lst[0]) and (atom_names[i] == "not_hydrogen")):
            listA.append(coordinates[i])
            residue_ids_A.append(residue_ids[i])
        elif ((chains[i] == lst[1]) and (atom_names[i] == "not_hydrogen")):
            listB.append(coordinates[i])
            residue_ids_B.append(residue_ids[i])      

    res_in_contact = []
    
    for i,posA in enumerate(listA):
        for j,posB in enumerate(listB):
            if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                cont = [residue_ids_A[i], residue_ids_B[j]]
                if cont not in res_in_contact:
                    res_in_contact.append(cont)
    
    return listA, listB, residue_ids_A, residue_ids_B, res_in_contact

def get_hydro_hits(protein):

    hydrophobicity_dict = {
            "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
            "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
            "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
            "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
        }

    hitsA = []
    hitsB = []

    listA, listB, resA, resB, res_list = external_contacts(protein)

    # print(f"Contacts in Chain A: {len(resA)} + Contacts in Chain B: {len(resB)}")
    # if len(resA) == len(resB):
    #     print("Lengths are equal. Proceeding to compare hydrophobicities.")
    # print("-----------------------------------------------------")

    for i in range(len(resA)):
        resnameA = get_residue_name(protein, resA[i])
    for i in range(len(resB)):
        resnameB = get_residue_name(protein, resB[i])
        
        if abs(hydrophobicity_dict[resnameA] - hydrophobicity_dict[resnameB]) <= 0.2:
            hitsA.append(resA[i])
            hitsB.append(resB[i])


    return hitsA, hitsB

def internal_contacts(protein):
    listA, listB, resA, resB, res_list = external_contacts(protein)
    
    lstA = {}
    lstB = {}
    
    for i,posA in enumerate(listA):
        tempA=[]
        for j,posB in enumerate(listA):
            if i != j:
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= 3.5:
                    tempA.append(resA[j])
        lstA[i] = tempA
    
    for i,posA in enumerate(listB):
        tempB=[]
        for j,posB in enumerate(listB):
            if i != j:
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= 3.5:
                    tempB.append(resB[j])
        lstB[i] = tempB

    return resA, resB, lstA, lstB
            
def create_graph(protein):
    hydroA = {}
    hydroB = {}
    
    resA, resB, graphA, graphB = internal_contacts(protein)
    listA, listB = get_hydro_hits(protein)

    for i in range(len(resA)):
        tempA =[]
        if resA[i] in listA:
            tempA.append(True)
        hydroA[i] = tempA

    for i in range(len(resB)):
        tempB =[]
        if resB[i] in listB:
            tempB.append(True)
        hydroB[i] = tempB

    return (graphA, graphB, hydroA, hydroB)
    # print(graphA)
    # print(graphB)
    # print(hydroA)
    # print(hydroB)

prot = Protein("/Users/smriti/Desktop/aeop/Protein-Decoy-Detection-SVR/targets/1c3a_complex_H.pdb")
#create_graph(prot)

lst = [(1, [1]), (2, [2, 1]), (3, [3, 2]), (4, [3, 2, 2]), (5, [4, 3, 5, 2])]

def score_fnc(lst):
    max_dist = len(lst)
    avg = 0
    score = 0
    SD = 0
    ns = []
    for dist in lst:
        dist_allowed, islands = dist
        tot = 0
        for island in islands:
            tot += len(islands) * island
        ns.append(tot)
        avg += tot
        score += tot/dist_allowed

    avg = avg/len(lst)
    for val in ns:
        SD += ((val - avg)**2)
    
    SD = SD/(max_dist-1)
    weighted_avg = score/max_dist

    return score, weighted_avg, avg, SD


'''
output_csv = f'/vast/palmer/scratch/ohern/sr2562/hydro_results/SD/{pdb_id}_hydrophobicity_fnc.csv'
        with open(output_csv, mode='w', newline='') as file:

            writer = csv.writer(file)
            writer.writerow(['pdb_file', 'avg_dif', 'N', 'standard_deviation'])
            for result in results:
                writer.writerow(result)
                
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
                prot = Protein(pdb_path)
                results.append((filename[:-4], (get_residues(prot))))
            else:
                print(f"File did not pass requirements.")

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
'''
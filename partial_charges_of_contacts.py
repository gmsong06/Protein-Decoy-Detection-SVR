from pdb_reader import Protein
import numpy as np
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
from pathlib import Path
import os
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

def get_residues(prot):
    coordinates = []
    atom_names = []                
    chains = []
    residue_ids = []
    atom_indices = []
    for model in prot.structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    chains.append(chain.id)
                    atom_names.append(atom.get_name())
                    coordinates.append(atom.coord)
                    residue_ids.append(residue.id[1])
                    atom_indices.append(atom.get_serial_number())

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

    listA = []
    listB = []
    residue_ids_A = []
    residue_ids_B = []
    atom_indices_A = []
    atom_indices_B = []
    lst = prot.get_chain_ids()
    for i in range(len(coordinates)):
        if ((chains[i] == lst[0]) and (atom_names[i] == "not_hydrogen")):
            listA.append(coordinates[i])
            residue_ids_A.append(residue_ids[i])
            atom_indices_A.append(atom_indices[i])
        elif ((chains[i] == lst[1]) and (atom_names[i] == "not_hydrogen")):
            listB.append(coordinates[i])
            residue_ids_B.append(residue_ids[i])
            atom_indices_B.append(atom_indices[i])      

    dist_thresh = 5
    res_in_contact = []
    for i, posA in enumerate(listA):
        for j, posB in enumerate(listB):
            if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                cont = [residue_ids_A[i], residue_ids_B[j], posA, posB, atom_indices_A[i], atom_indices_B[j]]
                if cont not in res_in_contact:
                    res_in_contact.append(cont)

    return res_in_contact

def calculate_atom_charges(pdb_file, atom_indices_in_contact):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    atom_charges = {atom_index: 0.0 for atom_index in atom_indices_in_contact}  # Initialize with zero charge

    # RDKit to compute MMFF94 partial charges
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    if mol:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200, nonBondedThresh=100.0)
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
        mmff_charges = mmff_props.GetMMFFPartialCharges()
        
        for atom in mol.GetAtoms():
            atom_index = atom.GetPDBResidueInfo().GetSerialNumber()
            if atom_index in atom_charges:
                atom_charges[atom_index] = mmff_charges[atom.GetIdx()]
    
    return atom_charges

def get_charge_similarity_score(file):
    protein = Protein(file)
    res_list = get_residues(protein)
    atom_indices_in_contact = {item[4] for item in res_list} | {item[5] for item in res_list}
    atom_charges = calculate_atom_charges(file, atom_indices_in_contact)
    similarity_score = 0

    print(len(res_list))
    for res in res_list:
        atom_indexA = res[4]
        atom_indexB = res[5]
        chargeA = atom_charges.get(atom_indexA, 0)
        chargeB = atom_charges.get(atom_indexB, 0)
        print(f"Contact between atoms {atom_indexA} and {atom_indexB}")
        print(f"Charges = {chargeA} and {chargeB}")
        
        # Score based on complementary charges
        similarity_score += 1 / (1 + np.abs(chargeA + chargeB))
        print(f"Similarity Score is {similarity_score}.")

    return similarity_score

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
                results.append((filename[:-4], get_charge_similarity_score(pdb_path)))
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

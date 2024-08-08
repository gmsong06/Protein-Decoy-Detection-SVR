import sys
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from scipy.spatial import KDTree
import csv
import argparse

def get_contacts(mol, cutoff=5.0):
    conf = mol.GetConformer()
    num_atoms = mol.GetNumAtoms()

    non_hydrogen_indices = []
    positions = []
    chains = []

    for i in range(num_atoms):
        atom_i = mol.GetAtomWithIdx(i)
        if atom_i.GetSymbol() not in {'H', 'h'}:
            non_hydrogen_indices.append(i)
            positions.append(conf.GetAtomPosition(i))
            chain_id = atom_i.GetPDBResidueInfo().GetChainId()
            chains.append(chain_id)

    positions = np.array(positions)
    tree = KDTree(positions)
    contacts = set()

    # Iterate over all non-hydrogen atoms
    for i, pos_i, chain_i in zip(non_hydrogen_indices, positions, chains):
        indices = tree.query_ball_point(pos_i, cutoff)
        for idx in indices:
            j = non_hydrogen_indices[idx]
            chain_j = chains[idx]
            if i < j and chain_i != chain_j:
                contacts.add((i, j))
    
    return list(contacts)

def test_mmff94_charges(pdb_file):
    try:
        mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if mol is None:
            raise ValueError("Molecule could not be read from the PDB file.")
    except Exception as e:
        return None

    mol = Chem.AddHs(mol)
    #print(f"Number of atoms after adding Hs: {mol.GetNumAtoms()}", flush=True)

    # Optimize the molecule
    optimize_status = AllChem.MMFFOptimizeMolecule(mol, maxIters=200, nonBondedThresh=100.0)


    try:
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94')
        if mmff_props is None:
            raise ValueError("MMFF properties could not be retrieved.")
    except Exception as e:
        #print(f"Error getting MMFF properties: {e}", flush=True)
        return None

    contacts = get_contacts(mol, cutoff=5.0)
    #print(f"Number of contacts: {len(contacts)}", flush=True)

    # Print partial charges of atoms in contact and calculate electric potential
    conf = mol.GetConformer()
    hits = 0 

    #print("Partial Charges of Atoms in Contact:", flush=True)
    for (i, j) in contacts:
        electric_potential = 0
        charge_i = mmff_props.GetMMFFPartialCharge(i)
        charge_j = mmff_props.GetMMFFPartialCharge(j)
        pos_i = conf.GetAtomPosition(i)
        pos_j = conf.GetAtomPosition(j)
        distance = pos_i.Distance(pos_j)
        electric_potential += (charge_i * charge_j) / distance
        if electric_potential < -0.05:
            hits +=1
            

    print(f"Total hits: {hits}", flush=True)
    return hits

def process_pdb_folder(full_folder_path, pdb_id):
    results = []
    relaxed_folder_path = os.path.join(full_folder_path, f"{pdb_id}_relaxed")
    random_folder_path = os.path.join(full_folder_path, f"random_negatives/rand_{pdb_id}_relaxed")

    paths = [relaxed_folder_path, random_folder_path]
    for path in paths:
        #print(f"Path is {path}")
        for filename in os.listdir(path):
            #print(f"Filename is {filename}")
            if filename.endswith('.pdb') and ("NoH" not in filename):
                pdb_path = os.path.join(path, filename)
                print(f"Processing {filename}", flush=True)
                score = test_mmff94_charges(pdb_path)
                if score is not None:
                    results.append((filename[:-4], score))
            else:
                print(f"File did not pass requirements.")

    output_csv = f'{pdb_id}_HITS.csv'
    with open(output_csv, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['pdb_file', 'partial_charges_score'])
        for result in results:
            writer.writerow(result)
            print(result)

def main(folder_path):
    for folder in os.listdir(folder_path):
        full_folder_path = os.path.join(folder_path, folder)
        if folder.startswith("sampled_") and os.path.isdir(full_folder_path):
            pdb_id = full_folder_path[-4:]
            #print(f"PDB id is {pdb_id}")
            process_pdb_folder(full_folder_path, pdb_id)
            #print("DONE----------------------------------------------------------------------")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
    args = parser.parse_args()
    main(args.pdb_folder)
import re
import Bio.PDB
from Bio.PDB.PDBIO import Select
import os
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from scipy.spatial import KDTree


class Protein:
    def __init__(self, pdb_file_path):
        self.pdb_file_path = pdb_file_path
        self.parser = Bio.PDB.PDBParser(QUIET=True)
        self.io = Bio.PDB.PDBIO()
        self.structure = self.parser.get_structure('protein', pdb_file_path)
        self.original_structure = self.structure.copy()  # Make a copy of the original structure
        self.hydrogen_removed = False

    def read_pdb(self):
        """
        Reads the PDB file and returns a list of lines
        
        """

        with open(self.pdb_file_path, 'r') as file:
            lines = file.readlines()
        return lines

    def get_chain_ids(self):
        """
        Returns two chain ids for the dimer structure

        """

        chain_id_1 = ''
        chain_id_2 = ''

        for model in self.structure:
            for chain in model:
                if not chain_id_1:
                    chain_id_1 = chain.id
                elif chain.id != chain_id_1:
                    chain_id_2 = chain.id
                    return chain_id_1, chain_id_2
        
        return [chain_id_1, chain_id_2]


    def get_atom_coords(self, chain_id: str):
        coords = []

        if chain_id in self.get_chain_ids():
            for residue in self.structure[0][chain_id]:
                for atom in residue:
                    coords.append({
                        atom.coord[0],
                        atom.coord[1],
                        atom.coord[2]
                    })
            return coords
        else:
            return "Invalid chain ID"
        
    
    def get_interface_residues(self):
        coordinates = []
        chains = []
        residue_ids = []

        # Collect data from the structure
        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    residue_id = residue.id[1]
                    for atom in residue:
                        chains.append(chain_id)
                        coordinates.append(atom.coord)
                        residue_ids.append(residue_id)

        # Separate coordinates and residue_ids by chain
        first_chain = chains[0]
        listA, listB = [], []
        residue_ids_A, residue_ids_B = [], []

        for i, chain in enumerate(chains):
            if chain == first_chain:
                listA.append(coordinates[i])
                residue_ids_A.append(residue_ids[i])
            else:
                listB.append(coordinates[i])
                residue_ids_B.append(residue_ids[i])

        # Find residues in contact
        dist_thresh = 5
        res_in_contact_A, res_in_contact_B = set(), set()

        for i, posA in enumerate(listA):
            for j, posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                    res_in_contact_A.add(residue_ids_A[i])
                    res_in_contact_B.add(residue_ids_B[j])

        # Combine and sort residues in contact
        all_residues_in_contact = sorted(res_in_contact_A | res_in_contact_B)

        return all_residues_in_contact

    def get_interface_atom_ids(self):
        # returns list of all the names of the ids of all atoms in the interface
        coordinates_dict = defaultdict(list)
        residue_ids_dict = defaultdict(list)
        chains = []

        # Collect data from the structure
        for model in self.structure:
            for chain in model:
                chain_id = chain.id
                for residue in chain:
                    for atom in residue:
                        chains.append(chain_id)
                        coordinates_dict[chain_id].append(atom.coord)
                        residue_ids_dict[chain_id].append(residue.id[1])

        # Get coordinates and residue IDs for the first chain and others
        first_chain = chains[0]
        listA = coordinates_dict[first_chain]
        residue_ids_A = residue_ids_dict[first_chain]

        atoms_in_contact = set()
        dist_thresh = 5

        for chain_id, listB in coordinates_dict.items():
            if chain_id == first_chain:
                continue
            
            residue_ids_B = residue_ids_dict[chain_id]

            for i, posA in enumerate(listA):
                for j, posB in enumerate(listB):
                    if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                        atoms_in_contact.add(i)
                        atoms_in_contact.add(j + len(residue_ids_A))

        # Convert the set to a sorted list
        return sorted(atoms_in_contact)

    
    def get_interface_atom_names(self):
        # returns names of atoms at interface
        coordinates = []
        atom_names = []
        chains = []
        residue_ids = []

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        chains.append(chain.id)
                        atom_names.append(atom.get_name())
                        coordinates.append(atom.coord)
                        residue_ids.append(residue.id[1])
 
        listA = []
        listB = []
        residue_ids_A = []
        residue_ids_B = []
        for i in range(len(coordinates)):
            if chains[i] == chains[0]:
                listA.append(coordinates[i])
                residue_ids_A.append(residue_ids[i])
            elif chains[i] != chains[0]:
                listB.append(coordinates[i])
                residue_ids_B.append(residue_ids[i])

        # finds atoms at interface and sorts them based on chain        
        atoms_in_contact = []
        A_positions = []
        B_positions = []
        for i,posA in enumerate(listA):
            for j,posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= 5:
                    A_position = i
                    A_positions.append(A_position)
                    B_position = j + len(residue_ids_A)
                    B_positions.append(B_position)
        A_positions_unique = list(set(A_positions))
        B_positions_unique = list(set(B_positions))
        for atom in A_positions_unique:
            atom_name = atom_names[atom]
            atoms_in_contact.append(atom_name)
        for atom in B_positions_unique:
            atom_name = atom_names[atom]
            atoms_in_contact.append(atom_name)

        return atoms_in_contact
    

    def get_interface_charges(self):
        # Initialize data structures
        coordinates = []
        chains = []
        amino_acid_names = []
        residue_ids = []
        residue_dict = {}

        # Collect data from the structure
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    residue_id = residue.id[1]
                    for atom in residue:
                        chains.append(chain.id)
                        coordinates.append(atom.coord)
                        amino_acid_names.append(residue.resname)
                        residue_ids.append(residue_id)
                        if residue_id not in residue_dict:
                            residue_dict[residue_id] = len(chains) - 1

        # Separate coordinates and residue_ids by chain
        first_chain = chains[0]
        listA, listB = [], []
        residue_ids_A, residue_ids_B = [], []

        for i, chain in enumerate(chains):
            if chain == first_chain:
                listA.append(coordinates[i])
                residue_ids_A.append(residue_ids[i])
            else:
                listB.append(coordinates[i])
                residue_ids_B.append(residue_ids[i])

        # Find residues in contact
        dist_thresh = 5
        res_in_contact_A, res_in_contact_B = set(), set()

        for i, posA in enumerate(listA):
            for j, posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                    res_in_contact_A.add(residue_ids_A[i])
                    res_in_contact_B.add(residue_ids_B[j])

        # Combine and sort residues in contact
        all_residues_in_contact = sorted(res_in_contact_A | res_in_contact_B)

        # Find indexes for the residues in contact
        indexes = [residue_dict[res_id] for res_id in all_residues_in_contact]

        # Define charges for amino acids
        charge_dict = {
            "LYS": +1, "ARG": +1, "HIS": 0.1,
            "ASP": -1, "GLU": -1
        }

        # Calculate charges in contact
        charges_in_contact = [
            charge_dict.get(amino_acid_names[index], 0) for index in indexes
        ]

        return charges_in_contact

            
    def get_hydrophobicities(self):
        # Initialize data structures
        coordinates = []
        chains = []
        amino_acid_names = []
        residue_ids = []
        residue_dict = {}

        # Collect data from the structure
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        chain_id = chain.id
                        residue_id = residue.id[1]
                        chains.append(chain_id)
                        coordinates.append(atom.coord)
                        amino_acid_names.append(residue.resname)
                        residue_ids.append(residue_id)
                        if residue_id not in residue_dict:
                            residue_dict[residue_id] = len(chains) - 1

        # Separate coordinates and residue_ids by chain
        first_chain = chains[0]
        listA, listB = [], []
        residue_ids_A, residue_ids_B = [], []

        for i, chain in enumerate(chains):
            if chain == first_chain:
                listA.append(coordinates[i])
                residue_ids_A.append(residue_ids[i])
            else:
                listB.append(coordinates[i])
                residue_ids_B.append(residue_ids[i])

        # Find residues in contact
        dist_thresh = 5
        res_in_contact_A, res_in_contact_B = set(), set()

        for i, posA in enumerate(listA):
            for j, posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                    res_in_contact_A.add(residue_ids_A[i])
                    res_in_contact_B.add(residue_ids_B[j])

        # Combine and sort residues in contact
        all_residues_in_contact = sorted(res_in_contact_A | res_in_contact_B)

        # Find indexes for the residues in contact
        indexes = [residue_dict[res_id] for res_id in all_residues_in_contact]

        # Hydrophobicity dictionary
        hydrophobicity_dict = {
            "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
            "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
            "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
            "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
        }

        # Calculate hydrophobicities
        hydrophobicities = [hydrophobicity_dict[amino_acid_names[index]] for index in indexes]

        return hydrophobicities


    def contact_number(self):
        coordinates = []
        atom_names = []
        chains = []
        residue_ids = []

        for model in self.structure:
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

    def atoms_in_contact_atom_ids(self):
        coordinates = []
        atom_ids = []
        atom_names = []
        chains = []

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        chains.append(chain.id)
                        atom_ids.append(atom.get_serial_number())  # Collect atom IDs
                        atom_names.append(atom.get_name())
                        coordinates.append(atom.coord)

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

        for i in range(len(coordinates)):
            if atom_names[i] == "not_hydrogen":
                chain_dict[chains[i]].append(coordinates[i])
                atom_dict[chains[i]].append(atom_ids[i])

        chain_ids = list(chain_dict.keys())
        if len(chain_ids) < 2:
            return 0

        listA = np.array(chain_dict[chain_ids[0]])
        listB = np.array(chain_dict[chain_ids[1]])

        atomsA = atom_dict[chain_ids[0]]
        atomsB = atom_dict[chain_ids[1]]

        # Using KDTree for efficient distance calculations
        treeA = KDTree(listA)
        dist_thresh = 5
        atom_ids_in_contact = set()

        for i, posB in enumerate(listB):
            indices = treeA.query_ball_point(posB, dist_thresh)
            for idx in indices:
                atom_ids_in_contact.add(atomsA[idx])
                atom_ids_in_contact.add(atomsB[i])

        sorted_atom_ids_in_contact = sorted(atom_ids_in_contact)

        return sorted_atom_ids_in_contact

if __name__ == "__main__":
    protein = Protein('targets/1acb_complex_H.pdb')
    #print(protein.get_chain_ids())
    #print(protein.get_atom_coords('E'))
    #print(protein.get_interface_residues())
    #print(protein.get_interface_atom_ids())
    #print(protein.get_interface_atom_names()) 
    #print(protein.get_interface_charges())
    #print(protein.get_hydrophobicities())
    #print(protein.get_residue_of_interface_atoms())
    print(protein.contact_number())
    print(protein.atoms_in_contact_atom_ids())

import re
import Bio.PDB
from Bio.PDB.PDBIO import Select
import os
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count


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
                            
        dist_thresh = 5
        res_in_contact = []
        for i,posA in enumerate(listA):
            for j,posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                    cont = [residue_ids_A[i], residue_ids_B[j]]
                    res_in_contact.append(cont)

        res_in_contact_A = []
        res_in_contact_B = []
        all_residues_in_contact = []
        
        for i in range(len(res_in_contact)):
            res_in_contact_A.append(res_in_contact[i][0])
            res_in_contact_B.append(res_in_contact[i][1])
        

        set_res_in_contact_A = list(set(res_in_contact_A))
        set_res_in_contact_B = list(set(res_in_contact_B))
        
        for i in range(len(set_res_in_contact_A)):
            all_residues_in_contact.append(set_res_in_contact_A[i])
        for i in range(len(set_res_in_contact_B)):
            all_residues_in_contact.append(set_res_in_contact_B[i])

        all_residues_in_contact.sort()            
        return all_residues_in_contact



    def get_interface_atom_names(self):
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

   

        dist_thresh = 5
        atoms_in_contact = []
        A_positions = []
        B_positions = []
        for i,posA in enumerate(listA):
            for j,posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
                    A_position = i
                    A_positions.append(A_position)
                    B_position = j + len(residue_ids_A)
                    B_positions.append(B_position)
        A_positions_unique = list(set(A_positions))
        B_positions_unique = list(set(B_positions))
        for atom in A_positions_unique:
            atoms_in_contact.append(atom)
        for atom in B_positions_unique:
            atoms_in_contact.append(atom)
        atoms_in_contact.sort()
        
        return atoms_in_contact
    
    def get_interface_atom_names(self):
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

        dist_thresh = 5
        atoms_in_contact = []
        A_positions = []
        B_positions = []
        for i,posA in enumerate(listA):
            for j,posB in enumerate(listB):
                if np.linalg.norm(np.array(posA) - np.array(posB)) <= dist_thresh:
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
        atoms_in_contact.sort()
        
        return atoms_in_contact

if __name__ == "__main__":
    protein = Protein('targets/1acb_complex_H.pdb')
    # print(protein.get_chain_ids())
    # print(protein.get_atom_coords('E'))
    #print(protein.get_interface_residues())
    #print(protein.get_interface_atom_ids())
    print(protein.get_interface_atom_names()) 
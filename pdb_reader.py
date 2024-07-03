import re
import Bio.PDB
from Bio.PDB.PDBIO import Select
import os
import argparse

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
        pass

    def get_interface_atoms(self):
        pass

if __name__ == "__main__":
    protein = Protein('targets/1acb_complex_H.pdb')
    # print(protein.get_chain_ids())
    print(protein.get_atom_coords('E'))
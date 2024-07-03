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
        """Reads a PDB file and returns a list of lines."""
        with open(self.pdb_file_path, 'r') as file:
            lines = file.readlines()
        return lines

from Bio.PDB import PDBParser
import numpy as np
import argparse
import csv
import os
from multiprocessing import Pool, cpu_count
from scipy.spatial import KDTree



from pdb_reader import Protein
Protein1 = Protein("targets/1acb_complex_H.pdb")
print(len(Protein1.get_residue_of_interface_atoms()))
print(len(Protein1.get_interface_atom_names()))
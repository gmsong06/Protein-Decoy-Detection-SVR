
from Bio.PDB import PDBParser
import numpy as np
import argparse
import csv
import os
from multiprocessing import Pool, cpu_count
from scipy.spatial import KDTree
from pdb_reader import Protein



from pdb_reader import Protein
Protein1 = Protein("targets/1acb_complex_H.pdb")


for item in Protein1.atoms_in_contact_atom_ids():
    print("hi")


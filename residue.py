import re
import Bio.PDB
from Bio.PDB.PDBIO import Select
import os
import argparse
import numpy as np
from multiprocessing import Pool, cpu_count
from collections import defaultdict
import pandas as pd
import freesasa

class Residue:
    def __init__(self, residue, chain):
        self.residue = residue
        self.chain = chain
        self.interactions = []
    
    def get_name(self):
        return self.residue.get_resname()

    def get_id(self):
        return self.residue.get_id()[1]
    
    def get_chain(self):
        return self.chain
    
    def add_interaction(self, res):
        self.interactions.append(res)
    
    def get_interactions(self):
        return self.interactions

    def get_label(self):
        return self.get_name() + "_" + str(self.get_id())
    
    def __eq__(self, other):
        return self.get_id() == other.get_id() and self.chain == other.chain

    def __hash__(self):
        return hash((self.get_id(), self.chain))
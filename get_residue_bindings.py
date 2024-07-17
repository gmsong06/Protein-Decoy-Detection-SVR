import os
from pdb_reader import Protein
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_folder", help="Folder path containing pdb files")
args = parser.parse_args()

def main():
    for file in os.listdir(args.input_folder):
        protein = Protein(os.path.join(args.input_folder, file))

        protein.get_interface_residue_names()
        
        break


if __name__=="__main__":
    main()
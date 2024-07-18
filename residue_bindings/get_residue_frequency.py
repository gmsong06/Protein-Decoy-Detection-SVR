import os
from pdb_reader import Protein
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument("input_folder", help="Folder path containing pdb files")
args = parser.parse_args()


def main():
    amino_acids = [
        'ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU',
        'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR',
    ]
    
    res_dict = {aa: 0 for aa in amino_acids}

    for file in os.listdir(args.input_folder):
        print(f"Processing file {file}")
        protein = Protein(os.path.join(args.input_folder, file))
        
        freq = protein.get_residue_frequency()

        for f in freq:
            if f in amino_acids:
                res_dict[f] += freq[f]
    
    print(res_dict)

    output_file = "residue_frequencies.pkl"
    with open(output_file, 'wb') as f:
        pickle.dump(res_dict, f)


if __name__=="__main__":
    main()

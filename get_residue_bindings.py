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
        'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR'
    ]
    
    res_dict = {aa: {aa_inner: 0 for aa_inner in amino_acids} for aa in amino_acids}

    # print(res_dict)
    for file in os.listdir(args.input_folder):
        print(f"Processing file {file}")
        protein = Protein(os.path.join(args.input_folder, file))

        res_name_A, res_name_B = protein.get_interface_residue_names()
        
        for i in range(len(res_name_A)):
            A = res_name_A[i]
            B = res_name_B[i]

            if A in amino_acids and B in amino_acids:
                res_dict[A][B] += 1
                res_dict[B][A] += 1
            else:
                print(f"EITHER {A} OR {B} IS NOT STANDARD")
    
    print(res_dict)
    
    output_file = "residue_contacts.pkl"
    with open(output_file, 'wb') as f:
        pickle.dump(res_dict, f)


if __name__=="__main__":
    main()
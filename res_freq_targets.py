import pickle
from pdb_reader import Protein, Complex
import math
import os
import pandas as pd

def calculate_rmsd(dict1, dict2):
    all_keys = set(dict1.keys()).union(dict2.keys())
    
    squared_diffs = []
    
    for key in all_keys:
        value1 = dict1.get(key, 0)
        value2 = dict2.get(key, 0)
        squared_diffs.append((value1 - value2) ** 2)
    
    mean_squared_diff = sum(squared_diffs) / len(squared_diffs)
    rmsd = math.sqrt(mean_squared_diff)
    
    return rmsd


# Define the function to process a single PDB file
def process_pdb(pdb_file_path, output_pkl_path):
    # Initialize an empty dictionary to store the data
    pdb_dict = {}

    # Extract the target name from the file path (assumes the file name contains the target name)
    target = pdb_file_path.split('/')[-1].split('_')[0]

    # Process the PDB file
    print(f"Processing target: {target}")
    complex = Complex(pdb_file_path)
    pdb_dict[target] = {'interactions': complex.get_interaction_dict()}

    # Print the dictionary (optional)
    print(pdb_dict)

    # Save the dictionary to a pickle file
    with open(output_pkl_path, 'wb') as f:
        pickle.dump(pdb_dict, f)

# Define the function to analyze interactions from a pickle file
def analyze_interactions(pkl_file_path):
    with open(pkl_file_path, 'rb') as file:
        pdb_dict = pickle.load(file)

    res_freq = {}
    res_dict = {}
    for target in pdb_dict:
        for interaction in pdb_dict[target]['interactions']:
            res_dict[interaction[0]] = pdb_dict[target]['interactions'][interaction][0]
            res_dict[interaction[1]] = pdb_dict[target]['interactions'][interaction][1]

            if interaction[0] not in res_freq:
                res_freq[interaction[0]] = 1
            else:
                res_freq[interaction[0]] += 1

            if interaction[1] not in res_freq:
                res_freq[interaction[1]] = 1
            else:
                res_freq[interaction[1]] += 1

    res_freq = dict(sorted(res_freq.items(), key=lambda item: item[1]))
    print(res_freq)
    return res_freq

# Example usage

path = '/home/as4643/palmer_scratch/Decoys/Supersampled_structures/sampled_1acb/1acb_relaxed/'
correct_dict = {38: 1, 57: 1, 58: 1, 92: 1, 97: 1, 100: 1, 140: 1, 299: 1, 300: 1, 301: 1, 
         145: 1, 185: 1, 186: 1, 187: 1, 209: 1, 271: 1, 216: 1, 222: 1, 223: 1, 
         283: 2, 285: 2, 56: 2, 141: 2, 304: 2, 302: 2, 289: 2, 147: 2, 274: 2, 
         171: 2, 287: 2, 190: 2, 37: 3, 39: 3, 40: 3, 168: 3, 189: 3, 210: 3, 
         191: 4, 211: 4, 213: 4, 55: 5, 144: 5, 275: 5, 212: 5, 276: 6, 277: 6, 
         282: 7, 214: 7, 280: 8, 188: 8, 278: 9, 281: 10, 279: 16}

df = pd.read_csv('final_data_groups.csv')

for file in os.listdir(path):

    pdb_file_path = f'/home/as4643/palmer_scratch/Decoys/Supersampled_structures/sampled_1acb/1acb_relaxed/{file}'  # Replace with your actual PDB file path
    output_pkl_path = f'{file[:-3]}.pkl'     # Replace with your desired output pickle file path

    process_pdb(pdb_file_path, output_pkl_path)

    # Analyze the interactions from the saved pickle file
    other_dict = analyze_interactions(output_pkl_path)

    rmsd = calculate_rmsd(correct_dict, other_dict)

    pdb_file = file[-4:]
    print()





# 0.982607368881035
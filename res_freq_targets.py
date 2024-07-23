import pandas as pd
import matplotlib.pyplot as plt
import pickle
from pdb_reader import Protein, Complex

# Load the CSV file into a DataFrame
targets = pd.read_csv('naomi_spearman_dict.csv')

gotten_dict = True

if not gotten_dict:
    # Initialize an empty dictionary to store the data
    targets_dict = {}

    # Populate the dictionary with 'pdb' and 'avg_spearman'
    for index, row in targets.iterrows():
        targets_dict[row['pdb']] = {'spearman': row['avg_spearman']}

    # for target in targets_dict:
    #     print(f"Processing target: {target}")
    #     protein = Protein(f'targets/{target}_complex_H.pdb')
    #     targets_dict[target]['fa'] = protein.get_fa()
    #     targets_dict[target]['contacts'] = protein.get_interface_residues()

    for target in targets_dict:
        print(f"Processing target: {target}")
        complex = Complex(f'targets/{target}_complex_H.pdb')
        targets_dict[target]['interactions'] = complex.get_interaction_dict()

    # Print the dictionary (optional)
    print(targets_dict)

    # Save the dictionary to a pickle file
    with open('targets_dict.pkl', 'wb') as f:
        pickle.dump(targets_dict, f)

else:
    with open('targets_dict.pkl', 'rb') as file:
        targets_dict = pickle.load(file)

    res_freq = {}
    res_dict = {}
    for target in targets_dict:
        for interaction in targets_dict[target]['interactions']:
            res_dict[interaction[0]] = targets_dict[target]['interactions'][interaction][0]
            res_dict[interaction[1]] = targets_dict[target]['interactions'][interaction][1]

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
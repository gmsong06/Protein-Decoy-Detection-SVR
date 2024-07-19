import pickle

# Path to the pickle file
input_file = "residue_frequencies.pkl"

# Load the dictionary from the pickle file
with open(input_file, 'rb') as f:
    res_dict = pickle.load(f)

def get_relative_freq(dict):
    total_contacts = 0
    for key, value in dict.items():
        # print(f"Amino acid is {key}")
        for aa in value:
            total_contacts += value[aa]
    
    for key, value in dict.items():
        for inner_aa in value:
            value[inner_aa] = round(value[inner_aa] / total_contacts, 5)

    return dict
    # print(total_contacts)

print(res_dict)
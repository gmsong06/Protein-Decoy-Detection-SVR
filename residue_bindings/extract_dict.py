import pickle

# Path to the pickle file
input_file = "residue_contacts.pkl"

# Load the dictionary from the pickle file
with open(input_file, 'rb') as f:
    res_dict = pickle.load(f)

def get_relative_freq(dict):
    for key, value in dict.items():
        print(f"Amino acid is {key}")
        total_contacts = 0
        for aa in value:
            total_contacts += value[aa]
        
        print(f"It has {total_contacts} total contacts")


        for inner_aa in value:
            value[inner_aa] = round(value[inner_aa] / total_contacts, 5)
    return dict
    # print(total_contacts)

relative_freq_dict = get_relative_freq(res_dict)
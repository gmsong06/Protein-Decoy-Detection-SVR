import pickle
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
import seaborn as sns

input_file = "residue_contacts.pkl"

# Load the dictionary from the pickle file
with open(input_file, 'rb') as f:
    res_dict = pickle.load(f)


interactions = {}

for aa in res_dict:
    if aa == 'UNK':
        continue
    
    for inner_aa in res_dict[aa]:
        if inner_aa == 'UNK':
            continue;
        interactions[(aa, inner_aa)] = res_dict[aa][inner_aa]

amino_acids = sorted(set([aa for pair in interactions.keys() for aa in pair]))

freq_matrix = np.zeros((len(amino_acids), len(amino_acids)))

# for (aa1, aa2), freq in freq_matrix.itemset():
#     i = amino_acids.index(aa1)
#     j = amino_acids.index(aa2)
#     freq_matrix[i, j] = freq
#     freq_matrix[j, i] = freq

for i in range(freq_matrix.shape[0]):
    for j in range(freq_matrix.shape[1]):
        freq_matrix[i][j] = interactions[(amino_acids[i], amino_acids[j])]

freq_df = pd.DataFrame(freq_matrix, index=amino_acids, columns=amino_acids)

scaler = MinMaxScaler()
normalized_freq_matrix = scaler.fit_transform(freq_df)
normalized_freq_df = pd.DataFrame(normalized_freq_matrix, index=amino_acids, columns=amino_acids)

plt.figure(figsize=(18, 12))
sns.heatmap(normalized_freq_df, annot=True, cmap="viridis")
plt.title("Frequency of Amino Acid Bindings")
plt.show()
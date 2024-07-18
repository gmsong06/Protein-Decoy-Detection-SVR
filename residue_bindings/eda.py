import pandas as pd
import pickle
import extract_dict
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dictionary from the pickle file
input_file = "residue_contacts.pkl"
with open(input_file, 'rb') as f:
    res_dict = pickle.load(f)

df = pd.DataFrame(res_dict)

df = df.drop(index='UNK', columns='UNK')

# plt.figure(figsize=(18, 12))
# sns.heatmap(df, annot=True, cmap="viridis")
# plt.title("Frequency of Amino Acid Bindings")
# plt.savefig('amino_acid_binding_frequencies.png')

# Binding probabilities
total_bindings = df.sum(axis=1)
binding_probabilities = df.div(total_bindings, axis=0)

binding_probabilities.to_csv('binding_probabilities.csv', index=True)

most_common_binding = df.idxmax(axis=1)

print(f"Most common binding {most_common_binding}")

plt.figure(figsize=(18, 12))
sns.heatmap(binding_probabilities, annot=True, cmap="viridis")
plt.title("Frequency of Amino Acid Bindings")
plt.show()
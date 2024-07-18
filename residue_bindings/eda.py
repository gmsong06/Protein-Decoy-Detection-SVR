import pandas as pd
import pickle
import extract_dict
import matplotlib.pyplot as plt
import seaborn as sns

# Load the dictionary from the pickle file
input_file = "residue_contacts.pkl"
with open(input_file, 'rb') as f:
    res_dict = pickle.load(f)

# Create a DataFrame from the dictionary
df = pd.DataFrame.from_dict(extract_dict.get_relative_freq(res_dict), orient='index')

# Define RSASA values
rsasa_values = {
    "ARG": 0.72002943, "ASP": 0.75367063, "GLU": 0.87591947, "LYS": 1, "ASN": 0.67819213,
    "GLN": 0.72278272, "PRO": 0.65123555, "HIS":  0.48907553, "SER": 0.52365422, "THR": 0.47798833,
    "GLY": 0.46477639, "TYR": 0.21646225, "ALA": 0.30953653, "CYS": 0, "MET": 0.18184843,
    "TRP":  0.14290738, "VAL": 0.10992156, "PHE": 0.0814021, "LEU": 0.10211201, "ILE": 0.06280283
}

# Create a DataFrame for RSASA values
rsasa_df = pd.DataFrame.from_dict(rsasa_values, orient='index', columns=['RSASA'])

# Melt the DataFrame to long format
df_melted = df.reset_index().melt(id_vars='index', var_name='Contact_AA', value_name='Frequency')

# Rename columns for clarity
df_melted.rename(columns={'index': 'AA'}, inplace=True)

# Merge RSASA values for both amino acids
df_melted = df_melted.merge(rsasa_df, left_on='AA', right_index=True)
df_melted = df_melted.merge(rsasa_df, left_on='Contact_AA', right_index=True, suffixes=('_AA', '_Contact_AA'))

# Calculate RSASA differences
df_melted['RSASA_Difference'] = abs(df_melted['RSASA_AA'] - df_melted['RSASA_Contact_AA'])

# Display the updated DataFrame with RSASA differences
print(df_melted.head())

# Plotting
plt.figure(figsize=(10, 6))
sns.scatterplot(x='RSASA_Difference', y='Frequency', data=df_melted, facecolors='none', edgecolor='black')
plt.title('RSASA Difference vs Binding Frequency')
plt.xlabel('RSASA Difference')
plt.ylabel('Binding Frequency')
plt.show()  # Ensure the plot is shown

# Save the plot
plt.savefig("plot.png")

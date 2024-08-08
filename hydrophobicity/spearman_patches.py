import pandas as pd
import matplotlib.pyplot as plt

# Load the CSV files
patches = pd.read_csv('hydrophobicity/patches_spearman.csv')
spearman_dict = pd.read_csv("naomi_spearman_dict.csv")
contacts = pd.read_csv("hydrophobicity/contacts_spearman.csv")

# Adjust patches and contacts DataFrames
patches['spearman_corr'] *= -1  # Multiply spearman_corr by -1 for patches
contacts['spearman_corr'] *= -1  # Multiply spearman_corr by -1 for contacts
patches['pdb'] = patches['filename'].str.extract(r'([0-9a-z]{4})')  # Extract pdb code from filename
contacts['pdb'] = contacts['pdb'].str.extract(r'([0-9a-z]{4})')  # Extract pdb code from filename

# Merge the DataFrames, keeping only matching rows
merged_df = pd.merge(spearman_dict, pd.merge(contacts[['pdb', 'spearman_corr']], 
                                             patches[['pdb', 'spearman_corr']], 
                                             on='pdb', suffixes=('_contacts', '_patches')), 
                     on='pdb', how='inner')

# Create lists for plotting
pdb = merged_df['pdb'].tolist()  # Use pdb codes for the x-axis labels
hydro_spearman = merged_df['spearman_corr_patches'].tolist()
contacts_spearman = merged_df['spearman_corr_contacts'].tolist()
avg_spearman = merged_df['avg_spearman'].tolist()

# Plot the data
fig, ax = plt.subplots(dpi=150, figsize=(10, 4))
ax.scatter(pdb, hydro_spearman, s=10, label='Hydrophobicity', color='red')
ax.scatter(pdb, contacts_spearman, s=10, label='Contacts', color='blue', marker='s')
ax.plot(pdb, avg_spearman, linewidth=1.5, label='Average', color='black')

ax.legend(fontsize=7)
ax.grid(axis='y')
ax.set_xlabel('PDB', fontsize=5)
ax.set_xticks(pdb)
ax.set_xticklabels(pdb, rotation=90, fontsize=5)  # Rotate labels for better readability
ax.tick_params(axis='both', which='major', labelsize=5)
ax.set_ylabel('Spearman', fontsize=7)
ax.set_title('Spearman Correlation by Scoring Function', fontsize=8)

plt.tight_layout()  # Adjust layout to prevent clipping of tick labels
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set font style and size similar to the reference code
rcParams['font.family'] = ['Computer Modern', 'serif']
rcParams['mathtext.fontset'] = 'cm'

# Load the CSV files
patches = pd.read_csv('hydrophobicity/patches_spearman.csv')
spearman_dict = pd.read_csv("naomi_spearman_dict.csv")
contacts = pd.read_csv("hydrophobicity/contacts_spearman.csv")
hydro = pd.read_csv("hydrophobicity/hydro_spearman.csv")  # Replace with the correct file for hydro

# Adjust patches, contacts, and hydro DataFrames
patches['spearman_corr'] *= -1  # Multiply spearman_corr by -1 for patches
contacts['spearman_corr'] *= -1  # Multiply spearman_corr by -1 for contacts
hydro['spearman_corr'] *= -1  # Multiply spearman_corr by -1 for hydro

# Extract pdb codes from filenames if needed
patches['pdb'] = patches['filename'].str.extract(r'([0-9a-z]{4})')  # Extract pdb code from filename
contacts['pdb'] = contacts['pdb'].str.extract(r'([0-9a-z]{4})')
hydro['pdb'] = hydro['pdb'].str.extract(r'([0-9a-z]{4})')  # Extract pdb code from pdb

# Merge the DataFrames, keeping only matching rows
merged_df = pd.merge(spearman_dict, 
                     pd.merge(contacts[['pdb', 'spearman_corr']], 
                              pd.merge(patches[['pdb', 'spearman_corr']], 
                                       hydro[['pdb', 'spearman_corr']], 
                                       on='pdb', suffixes=('_patches', '_hydro')),
                              on='pdb', suffixes=('_contacts', '_hydro_patches')),
                     on='pdb', how='inner')

# Create lists for plotting
pdb = merged_df['pdb'].tolist()  # Use pdb codes for the x-axis labels
patches_spearman = merged_df['spearman_corr_patches'].tolist()
contacts_spearman = merged_df['spearman_corr'].tolist()
avg_spearman = merged_df['avg_spearman'].tolist()

# Plot the data with the updated style
fig, ax = plt.subplots(dpi=150, figsize=(8, 4))  # Reduced figure size and DPI
ax.scatter(pdb, patches_spearman, s=50, label='Patches', marker='d', edgecolor='red', facecolors='none', linewidth=0.5)
ax.scatter(pdb, contacts_spearman, s=50, label='Contacts', marker='s', edgecolor='blue', facecolors='none', linewidth=0.5)
ax.plot(pdb, avg_spearman, linewidth=1.5, label='Average', color='black')  # Slightly thicker line

# Customize legend and grid
ax.legend(fontsize=8)  # Smaller legend font size
ax.grid(axis='y', linestyle='dotted', linewidth=0.5)

# Set axis labels and title with smaller font size
ax.set_xlabel('PDB', fontsize=10)
ax.set_ylabel(r'$\rho$', fontsize=10)
ax.set_title('Spearman Correlation by Scoring Function', fontsize=10)

# Customize x-ticks
ax.set_xticks(pdb)
ax.set_xticklabels(pdb, rotation=45, fontsize=8)  # Rotate labels 45 degrees
ax.tick_params(axis='both', which='major', labelsize=8)

plt.tight_layout()  # Adjust layout to prevent clipping of tick labels
plt.show()

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Load the CSV file
df = pd.read_csv('final_data_groups_hydro.csv')

# Extract columns hydro and contacts
hydro = df['hydro']
contacts = df['contacts']

# Calculate Pearson correlation coefficient
pearson_corr, _ = pearsonr(hydro, contacts)
print(f'Pearson correlation coefficient: {pearson_corr}')

# Create a scatter plot with custom point style
plt.figure(figsize=(10, 6))
sns.scatterplot(x=hydro, y=contacts, edgecolor='black', facecolor='none')

# Add labels and title with larger font size
plt.xlabel('Hydrophobicity', fontsize=36)
plt.ylabel('Contacts', fontsize=36)
plt.title('Hydrophobicity vs Contacts', fontsize=36)

# Increase tick label size
plt.xticks(fontsize=36)
plt.yticks(fontsize=36)

# Show the plot
plt.savefig('hydro_contacts.png')
plt.show()

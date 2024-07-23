import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# Load the CSV file
df = pd.read_csv('final_data_bsa_actual.csv')

# Extract columns hydro and contacts
hydro = df['bsa']
contacts = df['DockQ']

pearson_corr, _ = pearsonr(hydro, contacts)
print(f'Pearson correlation coefficient: {pearson_corr}')

# Create a scatter plot with custom point style
plt.figure(figsize=(10, 6))
sns.scatterplot(x=hydro, y=contacts, edgecolor='black', facecolor='none')

# Add labels and title
plt.xlabel('Hydrophobicity')
plt.ylabel('Contacts')
plt.title('Hydrophobicity vs Contacts')

# Show the plot
plt.show()

# Finds correlation between contacts and hydrophobic contacts

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

df = pd.read_csv('data/svr/final_data.csv')

hydro = df['hydro']
contacts = df['contacts']

# Calculate Pearson
pearson_corr, _ = pearsonr(hydro, contacts)
print(f'Pearson correlation coefficient: {pearson_corr}')

plt.figure(figsize=(10, 6))
sns.scatterplot(x=hydro, y=contacts, edgecolor='black', facecolor='none')

plt.xlabel('Hydrophobicity', fontsize=36)
plt.ylabel('Contacts', fontsize=36)
plt.title('Hydrophobicity vs Contacts', fontsize=36)

plt.xticks(fontsize=36)
plt.yticks(fontsize=36)

plt.savefig('hydro_contacts.png')
plt.show()

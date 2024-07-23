import pandas as pd
import numpy as np
#import os
#import argparse
import matplotlib.pyplot as plt
#from scipy.stats import pearsonr, spearmanr

hydro = pd.read_csv('hydro_spearman.csv')
spearman_dict = pd.read_csv("naomi_spearman_dict.csv")

merge = pd.merge(spearman_dict, hydro, on='pdb', how='outer')

pdb = [] #PDB
hydro_spearman = [] #HYDRO
avg_spearman = [] #AVG SPEARMAN


result = merge.sort_values(by=['avg_spearman'])
print(result)

for i in range(84):
    i+=1
    pdb.append(i)
for index, row in result.iterrows():
    hydro_spearman.append(float(row['spearman_correlation'])*(-1))
    avg_spearman.append(row['avg_spearman']) 

fig, ax = plt.subplots(dpi=150, figsize=(40,2))
ax.scatter(pdb, hydro_spearman, s=10, label='Hydrophbicity', color='red')
ax.plot(pdb, avg_spearman, linewidth=1.5, label='Average', color='black')
# ax.plot(pdb,hydro_spearman,color='red')
# ax.plot(pdb,avg_spearman,color='black')
ax.legend(fontsize = 7)
ax.grid(axis = 'y')
ax.set_xlabel('PDB', fontsize=7)
ax.set_xticks(pdb)
ax.tick_params(axis='both', which='major', labelsize=5)
ax.set_ylabel('Spearman', fontsize=7)
ax.set_title('Spearman Correlation by Scoring Function', fontsize=8)

#plt.savefig("/vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/spearman_plots/spearman_hydro.png")
plt.show()
plt.close()
import pandas as pd
import numpy as np
#import os
#import argparse
import matplotlib.pyplot as plt
#from scipy.stats import pearsonr, spearmanr

hydro = pd.read_csv('/vast/palmer/scratch/ohern/sr2562/hydro_spearman.csv')
spearman_dict = pd.read_csv('/vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/naomi_spearman_dict.csv')

result = pd.merge(spearman_dict, hydro, on='pdb', how='outer')

pdb = [] #PDB
hydro_spearman = [] #HYDRO
avg_spearman = [] #AVG SPEARMAN

for i in range(84):
    i+=1
    pdb.append(i)
for index, row in result.iterrows():
    hydro_spearman.append(float(row['spearman_correlation'])*(-1))
    avg_spearman.append(row['avg_spearman']) 

fig, ax = plt.subplots(dpi=150, figsize=(10,5))
ax.scatter(pdb, hydro_spearman, label='Hydrophbicity', color='red')
ax.scatter(pdb, avg_spearman, label='Average', marker='s', color='black')
ax.plot(pdb,hydro_spearman,color='red')
ax.plot(pdb,avg_spearman,color='black')
ax.legend()
ax.grid(axis = 'y')
ax.set_xlabel('PDB', fontsize=12)
ax.se_xticks(pdb)
ax.set_ylabel('Spearman', fontsize=12)
ax.set_title('Spearman Correlation by Scoring Function')

plt.savefig("/vast/palmer/scratch/ohern/sr2562/Protein-Decoy-Detection-SVR/spearman_plots/spearman_hydro.png")
plt.close()
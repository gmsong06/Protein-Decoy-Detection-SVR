import pandas as pd
import numpy as np
#import os
#import argparse
import matplotlib.pyplot as plt
#from scipy.stats import pearsonr, spearmanr

hydro = pd.read_csv('/Users/smriti/Desktop/aeop/Protein-Decoy-Detection-SVR/hydro_spearman.csv')
spearman_dict = pd.read_csv('/Users/smriti/Desktop/aeop/Protein-Decoy-Detection-SVR/naomi_spearman_dict.csv')

result = pd.merge(hydro, spearman_dict, on='pdb', how='outer')

x = [] #PDB
y = [] #HYDRO
avg = [] #AVG SPEARMAN

for i in range(84):
    i+=1
    x.append(i)
for index, row in result.iterrows():
    y.append(row['hydro'])
    avg.append(row['spearman_dict']) 

fig, ax = plt.subplots(dpi=150, figsize=(5, 5))
ax.scatter(x, y, label='Hydrophbicity', color='red')
ax.scatter(x, avg, label='Average', marker='s', color='black')
ax.plot(x,y,color='red')
ax.plot(x,avg,color='black')
ax.legend()
ax.grid(axis = 'y')
ax.set_xlabel('Spearman', fontsize=12)
ax.set_ylabel('PDB', fontsize=12)
ax.set_title('Spearman Correlation by Scoring Function')

plt.show()
plt.close()
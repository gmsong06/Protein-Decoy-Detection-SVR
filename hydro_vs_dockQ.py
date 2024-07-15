import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

dockq = pd.read_csv("final_data_groups.csv")

def main(folder_path):
    for file in os.listdir(args.predictions_path): 
        df = pd.read_csv(file)
        pdb_id = file[:4]
        dockq_scores = []
        hydro_scores = []
        results = []

        for index, row in dockq.iterrows():
            dockq_scores.append(row['DockQ'])
        for index, row in df.iterrows():
            hydro_scores.append(row['hydrophobicity_contacts'])
        
        spearman_corr, spearman_p = spearmanr(dockq_scores, hydro_scores)
        
        # Append results to the list
        results.append({'pdb_id': pdb_id, 'spearman_correlation': spearman_corr})

        fig, ax = plt.subplots(dpi=150, figsize=(5, 5))
        ax.scatter(hydro_scores, dockq_scores, edgecolor='black', color='none')

        ax.set_xlabel('Contacts of similar Hydrophbicity', fontsize=12)
        ax.set_ylabel('DockQ', fontsize=12)
        ax.set_title(f'Spearman Correlation: {spearman_corr:.2f}', fontsize=12)

        plot_filename = os.path.join("/Users/smriti/Desktop/aeop/Protein-Decoy-Detection-SVR/spearman_plots/hydro", f"{pdb_id}.png")
        plt.savefig(plot_filename)
        plt.close(fig)  # Close the figure to avoid memory issues

if __name__ == "__main__":
    main(args.pdb_folder)
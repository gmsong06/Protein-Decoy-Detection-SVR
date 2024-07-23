import pandas as pd
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

parser = argparse.ArgumentParser()
parser.add_argument("pdb_folder", type=str, help="Path to the folder containing PDB files")
args = parser.parse_args()

dockq = pd.read_csv("final_data_bsa_actual.csv")

def main(folder_path):
    results = []
    for file in os.listdir(folder_path): 
        file_path = os.path.join(folder_path, file)
        df = pd.read_csv(file_path)
        pdb_id = file[6:9]
        dockq_scores = []
        hydro_scores = []
        print(pdb_id)
        for index, row in dockq.iterrows():
            if row["pdb_id"] == pdb_id:
                dockq_scores.append(row['DockQ'])
        for index, row in df.iterrows():
            hydro_scores.append(row['interface_rsm'])
        
        print(dockq_scores)
        spearman_corr, spearman_p = spearmanr(dockq_scores, hydro_scores)
        
        # Append results to the list
        results.append({'pdb': pdb_id, 'spearman_correlation': spearman_corr})

        fig, ax = plt.subplots(dpi=150, figsize=(5, 5))
        ax.scatter(hydro_scores, dockq_scores, edgecolor='black', color='none')

        ax.set_xlabel('RSM', fontsize=12)
        ax.set_ylabel('DockQ', fontsize=12)
        ax.set_title(f'Spearman Correlation: {spearman_corr:.2f}', fontsize=12)

        plot_filename = os.path.join("spearman_plots/rsm", f"{pdb_id}.png")
        plt.savefig(plot_filename)
        plt.close(fig)  # Close the figure to avoid memory issues

    results_df = pd.DataFrame(results)
    results_df.to_csv('rsm_corr.csv', index=False)
    
if __name__ == "__main__":
    main(args.pdb_folder)
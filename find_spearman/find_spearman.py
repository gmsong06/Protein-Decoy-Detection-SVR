import argparse
import pandas as pd
from scipy.stats import spearmanr
import numpy as np
import matplotlib.pyplot as plt
import os

# Argument parser
parser = argparse.ArgumentParser()
parser.add_argument("predictions_path", type=str, help="Path to predictions")
parser.add_argument("output_path", type=str, help="Path to output figure")
args = parser.parse_args()

def main():
    # Read the data
    df = pd.read_csv(args.predictions_path)
    df['pdb_id'] = df['pdb_file'].str[-3:]
    dockq = df["actual_DockQ"]
    pred = df["prediction"]

    pdb_dict = {pdb_id: df[df['pdb_id'] == pdb_id] for pdb_id in df['pdb_id'].unique()}
    
    print("pdb_ids:", pdb_dict.keys())
    
    # Initialize a list to store results
    results = []

    for pdb_id, pdb_df in pdb_dict.items():
        print(f"PDB ID: {pdb_id}")
        
        dockq_scores = []
        pred_scores = []

        for index, row in pdb_df.iterrows():
            dockq_scores.append(row['actual_DockQ'])
            pred_scores.append(row['prediction'])
        
        spearman_corr, spearman_p = spearmanr(dockq_scores, pred_scores)
        
        # Append results to the list
        results.append({'pdb_id': pdb_id, 'spearman_correlation': spearman_corr})

        fig, ax = plt.subplots(dpi=150, figsize=(5, 5))
        ax.scatter(pred_scores, dockq_scores, edgecolor='black', color='none')

        ax.set_xlabel('Predicted DockQ', fontsize=12)
        ax.set_ylabel('DockQ', fontsize=12)
        ax.set_title(f'Spearman Correlation: {spearman_corr:.2f}', fontsize=12)

        fig.savefig(f'{args.output_path}/{pdb_id}.png')
        plt.close(fig)  # Close the figure to avoid memory issues

    # Convert results to DataFrame and save to CSV
    results_df = pd.DataFrame(results)
    results_df.to_csv(f'{args.output_path}/spearman_correlations.csv', index=False)

if __name__ == "__main__":
    main()

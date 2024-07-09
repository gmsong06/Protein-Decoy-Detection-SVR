import argparse
import pandas as pd
from scipy.stats import spearmanr

parser = argparse.ArgumentParser()
parser.add_argument("predictions_path", type=str, help="Path to predictions")
args = parser.parse_args()

def main():
    df = pd.read_csv(args.predictions_path)
    df['pdb_id'] = df['pdb_file'].str[-4:]

    pdb_dict = {pdb_id: df[df['pdb_id'] == pdb_id] for pdb_id in df['pdb_id'].unique()}
    
    print("pdb_ids:", pdb_dict.keys())
    
    for pdb_id, pdb_df in pdb_dict.items():
        print(f"PDB ID: {pdb_id}")
        
        dockq_scores = []
        pred_scores = []

        for index, row in pdb_df.iterrows():
            dockq_scores.append(row['actual_DockQ'])
            pred_scores.append(row['prediction'])

        
        spearman_corr, spearman_p = spearmanr(dockq_scores, pred_scores)


if __name__ == "__main__":
    main()

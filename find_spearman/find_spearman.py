import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("predictions_path", type=str, help="Path to predictions")
args = parser.parse_args()

def main():
    df = pd.read_csv(args.predictions_path)
    df['pdb_id'] = df['pdb_file'].str[-4:]

    pdb_dict = {pdb_id: df[df['pdb_id'] == pdb_id] for pdb_id in df['pdb_id'].unique()}
    
    print("Unique pdb_ids:", pdb_dict.keys())

    for pdb_id, pdb_df in pdb_dict.items():
        pdb_df.to_csv(f"{pdb_id}.csv", index=False)

if __name__ == "__main__":
    main()

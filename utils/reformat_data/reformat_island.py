import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("data_folder", type=str, help="Path to data folder")
args = parser.parse_args()

capri = True

def main():
    for data in os.listdir(args.data_folder):
        if data.endswith('.csv'):
            csv_path = os.path.join(args.data_folder, data)
            df = pd.read_csv(csv_path)

            df[['patch_alignment_score', 'chain_patch_consistency_score', 'avg_patch_score']] = df['patch_alignment_score'].str.extract(r'\(([^,]+), ([^,]+), ([^,]+)\)')

            # Convert the new columns to float
            df['patch_alignment_score'] = df['patch_alignment_score'].astype(float) * 100
            df['chain_patch_consistency_score'] = df['chain_patch_consistency_score'].astype(float) * 100
            df['avg_patch_score'] = df['avg_patch_score'].astype(float) * 100


            # Save the modified DataFrame to a new CSV file
            df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    main()
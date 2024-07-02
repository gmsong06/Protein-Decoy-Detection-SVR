import argparse
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("data_folder", type=str, help="Path to data folder")
args = parser.parse_args()

def main():
    for data_category in os.listdir(args.data_folder):
        category_path = os.path.join(args.data_folder, data_category)
        if os.path.isdir(category_path):
            for data in os.listdir(category_path):
                if data.endswith('.csv'):
                    csv_path = os.path.join(category_path, data)
                    pdb_id = data[:4]
                    print(f"PDB ID: {pdb_id}")
                    df = pd.read_csv(csv_path)

                    df['pdb_file'] = df['pdb_file'] + "_" + pdb_id

                    print(df['pdb_file'])

                    df.to_csv(csv_path, index=False)
if __name__ == "__main__":
    main()

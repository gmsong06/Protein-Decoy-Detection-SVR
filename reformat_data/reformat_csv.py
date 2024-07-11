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

            if capri:
                pdb_id = data[:3]
            else:
                pdb_id = data[:4]
            
            print(pdb_id)
            df = pd.read_csv(csv_path)

            print(df)
            for index, value in df['model'].items():
                if not value.endswith(pdb_id):
                    df.at[index, 'model'] = value + "_" + pdb_id
                
                if "random" in data:
                    if value.startswith("sampled"):
                        df.at[index, 'model'] = "random" + value[len("sampled"):].lstrip()
                    if value.startswith("relaxed"):
                        df.at[index, 'model'] = "random" + value[len("relaxed"):].lstrip()
                
                if "relaxed" in data:
                    if value.startswith("sampled"):
                        df.at[index, 'model'] = "relaxed" + value[len("sampled"):].lstrip()
                    if value.startswith("random"):
                        df.at[index, 'model'] = "relaxed" + value[len("random"):].lstrip()

            df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    main()
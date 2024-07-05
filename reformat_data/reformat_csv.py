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
                    df = pd.read_csv(csv_path)

                    for index, value in df['pdb_file'].items():
                        if value == 'random_3gfk_model_253_3gfk':
                            print("hi")
                        if not value.endswith(f"_corrected_H_0001_{pdb_id}"):
                            print(pdb_id)
                            df.at[index, 'pdb_file'] = value[:-5] + "_corrected_H_0001" + "_" + pdb_id
                        
                        if "random" in data:
                            if value.startswith("sampled"):
                                df.at[index, 'pdb_file'] = "random" + value[len("sampled"):].lstrip()
                            if value.startswith("relaxed"):
                                df.at[index, 'pdb_file'] = "random" + value[len("relaxed"):].lstrip()
                        
                        if "relaxed" in data:
                            if value.startswith("sampled"):
                                df.at[index, 'pdb_file'] = "relaxed" + value[len("sampled"):].lstrip()
                            if value.startswith("random"):
                                df.at[index, 'pdb_file'] = "relaxed" + value[len("random"):].lstrip()

    
                    df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    main()
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

            # if capri:
            #     pdb_id = data[len("capri_"): len("capri_") + 3]
            # else:
            pdb_id = data[:3]
            
            print(pdb_id)
            df = pd.read_csv(csv_path)

            for index, value in df['model'].items():
                if not value.endswith(f"_corrected_H_0001_{pdb_id}"):
                    print(pdb_id)
                    df.at[index, 'model'] = value[:-5] + "_corrected_H_0001" + "_" + pdb_id
                
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
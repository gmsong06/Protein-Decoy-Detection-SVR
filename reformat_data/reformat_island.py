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

            df[['score', 'weighted_avg', 'avg', 'standard_dev']] = df['score'].str.extract(r'\(([^,]+), ([^,]+), ([^,]+), ([^,]+)\)')

            # Convert the new columns to float
            df['score'] = df['score'].astype(float)
            df['weighted_avg'] = df['weighted_avg'].astype(float)
            df['avg'] = df['avg'].astype(float)
            df['standard_dev'] = df['standard_dev'].astype(float)


            # Save the modified DataFrame to a new CSV file
            df.to_csv(csv_path, index=False)


if __name__ == "__main__":
    main()
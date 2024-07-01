import os
import pandas as pd

def combine_data(dataframes):
    if dataframes:
        combined_df = pd.concat(dataframes, ignore_index=True).drop_duplicates(subset='pdb_file')
        combined_df.to_csv('full_data.csv', index=False)
        print("All CSV files combined and saved as 'full_data.csv'.")
    else:
        print("No CSV files found in the directory.")



def main(dir: str):
    dataframes = []

    for file in os.listdir(dir):
        file_path = os.path.join(dir, file)

        if os.path.isfile(file_path) and file.endswith('.csv'):
            print(f"Processing CSV file: {file}")
            df = pd.read_csv(file_path)
            dataframes.append(df)

    

if __name__ == "__main__":
    main("all_csvs")

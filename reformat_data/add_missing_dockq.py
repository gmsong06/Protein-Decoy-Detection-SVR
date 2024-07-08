import csv
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument("full_data", type=str, help="Path to full data")
parser.add_argument("missing_values", type=str, help="Path to additional DockQ scores")
args = parser.parse_args()


def fill_missing_values(file_path_full_data, file_path_missing_values):
    # Read the CSV files into DataFrames
    full_data = pd.read_csv(file_path_full_data)
    missing_values = pd.read_csv(file_path_missing_values)
    
    # Create a dictionary from missing_values for quick lookup
    missing_dict = dict(zip(missing_values['Decoy'], missing_values['DockQ']))
    
    # Loop through the rows in full_data and fill missing 'DockQ' values
    for index, row in full_data.iterrows():
        if pd.isna(row['DockQ']):
            pdb_file = row['pdb_file']
            if pdb_file in missing_dict:
                full_data.at[index, 'DockQ'] = missing_dict[pdb_file]

    return full_data


if __name__=="__main__":
    data = fill_missing_values(args.full_data, args.missing_values)
    data.to_csv('final_data.csv', index=False)
    print("Added to csv")

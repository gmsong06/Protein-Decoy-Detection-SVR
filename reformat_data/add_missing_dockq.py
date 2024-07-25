import csv
import argparse
import pandas as pd
import os


parser = argparse.ArgumentParser()
parser.add_argument("full_data", type=str, help="Path to full data")
parser.add_argument("missing_values", type=str, help="Path to folder containing additional DockQ scores")
args = parser.parse_args()


def fill_missing_values(file_path_full_data, folder_path_missing_values):
    # Read the full data CSV file into a DataFrame
    full_data = pd.read_csv(file_path_full_data)
    
    # Initialize an empty dictionary to store missing values
    missing_dict = {}

    # Loop through all files in the missing_values folder
    for file_name in os.listdir(folder_path_missing_values):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path_missing_values, file_name)
            missing_values = pd.read_csv(file_path)
            
            # Update the dictionary with values from the current file
            missing_dict.update(dict(zip(missing_values['Decoy'], missing_values['DockQ'])))

    # Loop through the rows in full_data and fill missing 'DockQ' values
    for index, row in full_data.iterrows():
        print(row)
        if pd.isna(row['DockQ']):
            pdb_file = row['pdb_file']

            if pdb_file in missing_dict:
                # print(missing_dict)
                full_data.at[index, 'DockQ'] = missing_dict[pdb_file]

    return full_data


if __name__ == "__main__":
    data = fill_missing_values(args.full_data, args.missing_values)
    data.to_csv('final_data_capri_with_hydro.csv', index=False)
    print("Added to csv")

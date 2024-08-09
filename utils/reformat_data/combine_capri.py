# Not in use

import pandas as pd
import os


def main(folder_path):
    # Get all files in the folder
    files = os.listdir(folder_path)

    # Dictionary to store pairs of CSV files by PDB ID
    pdb_files = {}

    for file in files:
        if file.endswith('_random_dockQ.csv') or file.endswith('_sampled_dockQ.csv'):
            pdb_id = file.split('_')[0]
            if pdb_id not in pdb_files:
                pdb_files[pdb_id] = {'random': None, 'sampled': None}
            if file.endswith('_random_dockQ.csv'):
                pdb_files[pdb_id]['random'] = file
            elif file.endswith('_sampled_dockQ.csv'):
                pdb_files[pdb_id]['sampled'] = file

    # Process each PDB ID
    for pdb_id, csv_files in pdb_files.items():
        random_file = csv_files['random']
        sampled_file = csv_files['sampled']

        if random_file and sampled_file:
            # Read the CSV files
            random_df = pd.read_csv(os.path.join(folder_path, random_file))
            sampled_df = pd.read_csv(os.path.join(folder_path, sampled_file))

            # Combine the DataFrames
            combined_df = pd.concat([random_df, sampled_df], ignore_index=True)

            # Save the combined DataFrame to a new CSV file
            combined_file_path = os.path.join(folder_path, f'{pdb_id}_combined_dockQ.csv')
            combined_df.to_csv(combined_file_path, index=False)
            print(f'Combined CSV for {pdb_id} saved as {combined_file_path}')
        else:
            print(f'Missing files for PDB ID {pdb_id}. Ensure both random and sampled files are present.')


if __name__ == "__main__":
    main("/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/capri_csvs/dockQ")

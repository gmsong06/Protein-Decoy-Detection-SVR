import os
import requests
import pandas as pd

# Path to the CSV file containing the list of PDB IDs and assembly IDs
pdb_list_csv = 'data/supersampled_candidates_list.csv'

# Directory to save the downloaded PDB files
output_dir = 'pdb_assembly_files'

# Create the output directory if it does not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Read the list of PDB IDs and assembly IDs from the CSV file
pdb_assembly_df = pd.read_csv(pdb_list_csv)

# Base URL for downloading PDB files for specific assemblies
base_url = 'https://files.rcsb.org/download/'

for index, row in pdb_assembly_df.iterrows():
    pdb_id = row['pdb_id']
    assembly_id = row['biounit_no']
    url = f'{base_url}{pdb_id}.pdb{assembly_id}.gz'
    response = requests.get(url)

    if response.status_code == 200:
        pdb_file_path = os.path.join(output_dir, f'{pdb_id}_assembly{assembly_id}.pdb.gz')
        with open(pdb_file_path, 'wb') as pdb_file:
            pdb_file.write(response.content)
        print(f'Downloaded {pdb_id}_assembly{assembly_id}.pdb.gz')
    else:
        print(f'Failed to download {pdb_id}_assembly{assembly_id}.pdb.gz')

print('Download completed.')

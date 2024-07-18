import pandas as pd

# Load the data from CSV files
dockq_scores = pd.read_csv('final_data_groups.csv')
to_be_added = pd.read_csv('final_data_hydro_groups.csv')

# Merge the DataFrames on the 'pdb_file' column
merged_df = pd.merge(to_be_added, dockq_scores[['pdb_file', 'DockQ']], on='pdb_file', how='left')

# Save the merged DataFrame back to a CSV file (optional)
merged_df.to_csv('merged_data.csv', index=False)

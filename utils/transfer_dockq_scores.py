# Not in use

import pandas as pd

dockq_scores = pd.read_csv('final_data_capri_groups.csv')
to_be_added = pd.read_csv('final_data_bsa.csv')

merged_df = pd.merge(to_be_added, dockq_scores[['pdb_file', 'DockQ']], on='pdb_file', how='left')

merged_df.to_csv('merged_data.csv', index=False)

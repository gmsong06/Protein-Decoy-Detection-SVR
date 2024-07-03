import pandas as pd
import csv
import os


main_file = pd.read_csv("combined_data.csv", on_bad_lines='skip')

def dockq(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith('.csv') and os.path.isdir(file_path):
            data = pd.read_csv(file_path, on_bad_lines='skip')
            pdb_id = filename[:4]

            for index, value in data['Decoys'].items():
                        if data.at[index, 'Decoys'][-4:] != pdb_id:
                            data.at[index, 'Decoys'] = data.at[index, 'Decoys'] + "_" + pdb_id
                            print("Appending " + pdb_id)

            data.to_csv(file_path, index=False)
            pd.concat([main_file, data["DockQ"]],
                  axis = 1)
                        
    

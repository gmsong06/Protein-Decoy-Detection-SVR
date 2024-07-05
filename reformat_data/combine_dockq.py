import pandas as pd
import csv
import os


main_file = pd.read_csv("combined_data.csv", on_bad_lines='skip')
lst = {}
def dockQ(folder_path):
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if filename.endswith('.csv'):
            data = pd.read_csv(file_path, on_bad_lines='skip')
            pdb_id = filename[:4]

            for index, value in main_file['pdb_file'].items():
                  lst.add(value[-4:])

            for index, value in data['Decoy'].items():
                        if (not data.at[index, 'Decoy'].endswith(f"_corrected_H_0001_{pdb_id}")) and (pdb_id in lst):
                            data.at[index, 'Decoy'] = data.at[index, 'Decoy'] + "_corrected_H_0001" + "_" + pdb_id
                            print("Appending " + pdb_id)

            data.to_csv(file_path, index=False)
            pd.concat([main_file, data["DockQ"]],
                  axis = 1)
            main_file.to_csv('combined_data_2.csv', index=False)
                        

if __name__ == "__main__":
      dockQ('final_manuscript_supersample_scores')

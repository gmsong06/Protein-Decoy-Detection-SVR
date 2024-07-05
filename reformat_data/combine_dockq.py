import pandas as pd
import os

# Load the main file
main_file = pd.read_csv("combined_data.csv")

# Create a set of unique pdb_file suffixes
lst = set()
for value in main_file['pdb_file']:
    lst.add(value[-4:])

# Define the function to process DockQ scores
def dockQ(folder_path):
    count = 0
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        
        if filename.endswith('.csv') and filename[:4] in lst:
            data = pd.read_csv(file_path)
            pdb_id = filename[:4]

            # Append the suffix to the 'Decoy' column values if necessary
            for index, value in data['Decoy'].items():
                if not value.endswith(f"_corrected_H_0001_{pdb_id}"):
                    data.at[index, 'Decoy'] = value + "_corrected_H_0001" + "_" + pdb_id
                    print(f"Appending {pdb_id} to Decoy value at index {index}")

            print(f"PDB ID is {pdb_id}")

            # Save the modified file
            data.to_csv(file_path, index=False)

            # Align DockQ scores based on matching Decoy with pdb_file
            for index, row in main_file.iterrows():
                pdb_file = row['pdb_file']
                matching_rows = data[data['Decoy'] == pdb_file]
                if not matching_rows.empty:
                    
                    dockq_value = matching_rows['DockQ'].values[0]
                    print(f"DockQ value is {dockq_value}")
                    print(f"File is {row}")
                    main_file.at[index, 'DockQ'] = dockq_value

    # Save the combined data to a new file
    main_file.to_csv('combined_data_3.csv', index=False)

# Run the dockQ function
if __name__ == "__main__":
    dockQ('final_manuscript_supersample_scores')

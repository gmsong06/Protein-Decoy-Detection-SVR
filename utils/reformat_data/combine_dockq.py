# Not in use

import pandas as pd
import os

# Load the main file
main_file = pd.read_csv("combined_data.csv")

# Create a set of unique pdb_file prefixes
unique_pdb_prefixes = set(main_file['pdb_file'].str[-4:])


# Function to process DockQ scores
def process_dockq_scores(folder_path):
    combined_decoy_map = {}

    for filename in os.listdir(folder_path):
        if filename.endswith('.csv') and filename[:4] in unique_pdb_prefixes:
            file_path = os.path.join(folder_path, filename)
            data = pd.read_csv(file_path)
            pdb_id = filename[:4]

            # Append the suffix to the 'Decoy' column values and modify based on filename conditions
            def modify_decoy(decoy):
                if not decoy.endswith(f"_corrected_H_0001_{pdb_id}"):
                    decoy += f"_corrected_H_0001_{pdb_id}"

                if "random" in filename:
                    if decoy.startswith('sampled'):
                        decoy = "random" + decoy[len("sampled"):].lstrip()
                    elif decoy.startswith('relaxed'):
                        decoy = "random" + decoy[len("relaxed"):].lstrip()

                elif "relaxed" in filename:
                    if decoy.startswith('sampled'):
                        decoy = "relaxed" + decoy[len("sampled"):].lstrip()
                    elif decoy.startswith('random'):
                        decoy = "relaxed" + decoy[len("random"):].lstrip()

                return decoy

            data['Decoy'] = data['Decoy'].apply(modify_decoy)

            # Save the modified file
            data.to_csv(file_path, index=False)

            # Create a dictionary mapping Decoy to DockQ
            decoy_map = data.set_index('Decoy')['DockQ'].to_dict()

            # Update the combined decoy_map
            combined_decoy_map.update(decoy_map)

    # Update DockQ values in main_file based on combined_decoy_map
    main_file['DockQ'] = main_file['pdb_file'].apply(lambda x: combined_decoy_map.get(x, float('nan')))

    # Save the combined data to a new file
    main_file.to_csv('combined_data_3.csv', index=False)


# Run the process_dockq_scores function
if __name__ == "__main__":
    process_dockq_scores('final_manuscript_supersample_scores')

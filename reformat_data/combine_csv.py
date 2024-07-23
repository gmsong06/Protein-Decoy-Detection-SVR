import pandas as pd
import os

def main(full_dir: str):
    all_contacts_dfs = []
    interface_rsm_dfs = []
    interface_flatness_dfs = []
    hydrophobicity_dfs = []
    bsa_dfs = []

    # Traverse the directory and read the CSV files
    for data_category in os.listdir(full_dir):
        data_category_path = os.path.join(full_dir, data_category)
        if os.path.isdir(data_category_path):
            for data_csv in os.listdir(data_category_path):
                file_path = os.path.join(data_category_path, data_csv)
                if data_csv.endswith('_all_contacts.csv'):
                    print(f"Reading {file_path}")
                    all_contacts_dfs.append(pd.read_csv(file_path))
                elif data_csv.endswith('_interface_rsm.csv'):
                    print(f"Reading {file_path}")
                    interface_rsm_dfs.append(pd.read_csv(file_path))
                elif data_csv.endswith('_interface_flatness.csv'):
                    print(f"Reading {file_path}")
                    interface_flatness_dfs.append(pd.read_csv(file_path))
                elif data_csv.endswith('_hydrophobicity.csv'):
                    print(f"Reading {file_path}")
                    hydrophobicity_dfs.append(pd.read_csv(file_path))
                elif data_csv.endswith('_bsa.csv'):
                    print(f"Reading {file_path}")
                    bsa_dfs.append(pd.read_csv(file_path))

    # Check the lengths of the lists
    print(f"Found {len(all_contacts_dfs)} '_all_contacts.csv' files")
    print(f"Found {len(interface_rsm_dfs)} '_interface_rsm.csv' files")
    print(f"Found {len(interface_flatness_dfs)} '_interface_flatness.csv' files")
    print(f"Found {len(hydrophobicity_dfs)} '_hydrophobicity.csv' files")
    print(f"Found {len(bsa_dfs)} '_bsa.csv' files")

    combined_df = pd.DataFrame()

    # Merge the DataFrames
    for idx, all_contacts_df in enumerate(all_contacts_dfs):
        pdb_file = all_contacts_df['pdb_file'].iloc[0]
        interface_rsm_df = next((df for df in interface_rsm_dfs if pdb_file in df['pdb_file'].values), None)
        interface_flatness_df = next((df for df in interface_flatness_dfs if pdb_file in df['pdb_file'].values), None)
        # hydrophobicity_df = next((df for df in hydrophobicity_dfs if pdb_file in df['pdb_file'].values), None)
        bsa_df = next((df for df in bsa_dfs if pdb_file in df['pdb_file'].values), None)

        # Debugging information
        print(f"Processing set {idx+1}:")
        print(f"  PDB file: {pdb_file}")
        print(f"  Found interface_rsm_df: {'Yes' if interface_rsm_df is not None else 'No'}")
        print(f"  Found interface_flatness_df: {'Yes' if interface_flatness_df is not None else 'No'}")
        # print(f"  Found hydrophobicity_df: {'Yes' if hydrophobicity_df is not None else 'No'}")

        if interface_rsm_df is not None and interface_flatness_df is not None and bsa_df is not None:
            print(f"Merging set {idx+1}")
            try:
                merged_df = pd.merge(all_contacts_df, pd.merge(interface_rsm_df, pd.merge(interface_flatness_df, bsa_df, on='pdb_file'), on='pdb_file'), on='pdb_file')
                combined_df = pd.concat([combined_df, merged_df], ignore_index=True)
            except Exception as e:
                print(f"Error merging set {idx+1}: {e}")
        else:
            print(f"Skipping set {idx+1} due to missing files")

    # Check if combined_df is empty before saving
    if not combined_df.empty:
        combined_df.to_csv('combined_data.csv', index=False)
        print("CSV files have been successfully combined into 'combined_data.csv'")
    else:
        print("No data to combine. Please check the input files.")

if __name__ == "__main__":
    main('capri_csvs')

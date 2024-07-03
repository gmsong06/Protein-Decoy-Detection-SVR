import pandas as pd
import os

def main(full_dir: str):
    all_contacts_dfs = []
    interface_rsm_dfs = []
    
    for data_category in os.listdir(full_dir):
        data_category_path = os.path.join(full_dir, data_category)
        if os.path.isdir(data_category_path):
            for data_csv in os.listdir(data_category_path):
                if data_csv.endswith('_all_contacts.csv'):
                    all_contacts_dfs.append(pd.read_csv(os.path.join(data_category_path, data_csv)))
                elif data_csv.endswith('_interface_rsm.csv'):
                    interface_rsm_dfs.append(pd.read_csv(os.path.join(data_category_path, data_csv)))
    
    combined_df = pd.DataFrame()
    
    for all_contacts_df, interface_rsm_df in zip(all_contacts_dfs, interface_rsm_dfs):
        merged_df = pd.merge(all_contacts_df, interface_rsm_df, on='pdb_file')
        combined_df = pd.concat([combined_df, merged_df], ignore_index=True)
    
    combined_df.to_csv('combined_data.csv', index=False)
    print("CSV files have been successfully combined into 'combined_interface_rsm.csv'")

if __name__ == "__main__":
    main('all_csvs')

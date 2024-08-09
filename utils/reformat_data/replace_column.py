# Not in use

import os
import pandas as pd


def merge_csv_files(directory, output_file):
    all_dfs = []

    # Traverse the directory and read the CSV files
    for file_name in os.listdir(directory):
        file_path = os.path.join(directory, file_name)
        if file_name.endswith('.csv'):
            print(f"Reading {file_path}")
            try:
                df = pd.read_csv(file_path)
                all_dfs.append(df)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")

    # Concatenate all DataFrames
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)
        try:
            combined_df.to_csv(output_file, index=False)
            print(f"All CSV files have been successfully combined into '{output_file}'")
        except Exception as e:
            print(f"Error saving combined data to CSV: {e}")
    else:
        print("No CSV files found or all CSV files are empty")


def merge_with_additional_csv(combined_csv_file, additional_csv_file, output_file):
    try:
        combined_df = pd.read_csv(combined_csv_file)
        additional_df = pd.read_csv(additional_csv_file)
    except Exception as e:
        print(f"Error reading CSV files: {e}")
        return

    # Merge the DataFrames on the 'pdb_file' column
    try:
        merged_df = pd.merge(combined_df, additional_df, on='pdb_file', how='left')
        merged_df.to_csv(output_file, index=False)
        print(f"Merged data has been saved to '{output_file}'")
    except Exception as e:
        print(f"Error merging DataFrames: {e}")


if __name__ == "__main__":
    directory = 'hydro_norm_temp'  # Replace with the path to your directory
    combined_output_file = 'combined_norm_hydro.csv'  # Replace with your desired output file name

    additional_csv_file = 'final_data_capri_with_hydro.csv'  # Replace with the path to your additional CSV file
    final_output_file = 'final_merged_data.csv'  # Replace with your desired final output file name

    # Merge all CSV files in the directory
    merge_csv_files(directory, combined_output_file)

    # Merge the combined CSV with the additional CSV
    merge_with_additional_csv(combined_output_file, additional_csv_file, final_output_file)

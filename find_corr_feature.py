import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

# Load the full_data DataFrame
full_data = pd.read_csv('final_data_groups_hydro.csv')

def process_folder(folder_path, full_data, output_folder, spearman_output_file):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # List to store Spearman correlation results
    spearman_results = []

    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            print(filename)
            file_path = os.path.join(folder_path, filename)
            # Load the current CSV file
            data = pd.read_csv(file_path)
            
            dockQ_values = []
            weighted_avg_values = []

            for _, row in data.iterrows():
                pdb_file = row['pdb_file']
                weighted_avg = row['avg']
                
                # Find the matching row in full_data
                matching_row = full_data[full_data['pdb_file'] == pdb_file]
                
                if not matching_row.empty:
                    dockQ = matching_row['DockQ'].values[0]
                    dockQ_values.append(dockQ)
                    weighted_avg_values.append(weighted_avg)
    
            # Calculate Spearman correlation
            if dockQ_values and weighted_avg_values:
                spearman_corr, _ = spearmanr(dockQ_values, weighted_avg_values)
                spearman_results.append({'filename': filename, 'spearman_corr': spearman_corr})
            
            # Plotting dockQ vs weighted_avg for the current file
            print("Plotting")
            plt.figure(figsize=(10, 6))
            plt.scatter(dockQ_values, weighted_avg_values, facecolors='none', edgecolors='black')
            plt.title(f'DockQ vs weighted_avg for {filename}')
            plt.xlabel('DockQ')
            plt.ylabel('weighted_avg')
            plt.grid(False)
            
            # Save the plot
            plot_filename = os.path.join(output_folder, f'{filename}.png')
            plt.savefig(plot_filename)
            print("Saved plot")
            plt.close()
            print()
    
    # Save Spearman results to CSV
    spearman_df = pd.DataFrame(spearman_results)
    spearman_df.to_csv(spearman_output_file, index=False)
    print(f"Spearman correlations saved to {spearman_output_file}")

# Specify the folder containing the CSV files, the output folder for plots, and the output file for Spearman correlations
csv_folder_path = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/all_csvs/islands'
output_folder_path = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/spearman_plots/islands/avg'
spearman_output_file = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/spearman_plots/islands/avg/spearman_results.csv'

# Process the folder, create plots, and save Spearman correlations
process_folder(csv_folder_path, full_data, output_folder_path, spearman_output_file)

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
                weighted_avg = row['avg_patch_score']
                
                # Find the matching row in full_data
                matching_row = full_data[full_data['pdb_file'] == pdb_file]
                
                if not matching_row.empty:
                    dockQ = matching_row['DockQ'].values[0]
                    dockQ_values.append(dockQ)
                    weighted_avg_values.append(weighted_avg)
    
            # Calculate Spearman correlation
            if dockQ_values and weighted_avg_values:
                spearman_corr, _ = spearmanr(weighted_avg_values, dockQ_values)
                spearman_results.append({'filename': filename, 'spearman_corr': spearman_corr})
            
            # Plotting avg_patch_score vs dockQ for the current file
            print("Plotting")
            plt.figure(figsize=(10, 6))
            plt.scatter(weighted_avg_values, dockQ_values, facecolors='none', edgecolors='black')
            plt.title(f'avg_patch_score vs DockQ for {filename}')
            plt.xlabel('avg_patch_score')
            plt.ylabel('DockQ')
            plt.grid(False)

            # Setting axis limits based on data range with a margin
            if weighted_avg_values:
                x_min, x_max = min(weighted_avg_values), max(weighted_avg_values)
                plt.xlim(x_min - 0.1*(x_max - x_min), x_max + 0.1*(x_max - x_min))
            if dockQ_values:
                y_min, y_max = min(dockQ_values), max(dockQ_values)
                plt.ylim(y_min - 0.1*(y_max - y_min), y_max + 0.1*(y_max - y_min))
            
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
output_folder_path = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/spearman_plots/islands/avg_patch_score'
spearman_output_file = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/spearman_plots/islands/avg_patch_score/spearman_results.csv'

# Process the folder, create plots, and save Spearman correlations
process_folder(csv_folder_path, full_data, output_folder_path, spearman_output_file)

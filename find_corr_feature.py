import os
import pandas as pd
import matplotlib.pyplot as plt

# Load the full_data DataFrame
full_data = pd.read_csv('final_data_groups_hydro.csv')

def process_folder(folder_path, full_data, output_folder):
    # Ensure the output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    # Loop through all files in the folder
    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            file_path = os.path.join(folder_path, filename)
            # Load the current CSV file
            data = pd.read_csv(file_path)
            
            dockQ_values = []
            weighted_avg_values = []

            for _, row in data.iterrows():
                pdb_file = row['pdb_file']
                weighted_avg = row['weighted_avg']
                
                # Find the matching row in full_data
                matching_row = full_data[full_data['pdb_file'] == pdb_file]
                
                if not matching_row.empty:
                    dockQ = matching_row['dockQ'].values[0]
                    dockQ_values.append(dockQ)
                    weighted_avg_values.append(weighted_avg)
    
            # Plotting dockQ vs weighted_avg for the current file
            plt.figure(figsize=(10, 6))
            plt.scatter(dockQ_values, weighted_avg_values, alpha=0.5)
            plt.title(f'DockQ vs weighted_avg for {filename}')
            plt.xlabel('dockQ')
            plt.ylabel('weighted_avg')
            plt.grid(True)
            
            # Save the plot
            plot_filename = os.path.join(output_folder, f'plot_{filename}.png')
            plt.savefig(plot_filename)
            plt.close()

# Specify the folder containing the CSV files and the output folder for plots
csv_folder_path = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/all_csvs/islands'
output_folder_path = '/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/spearman_plots/islands'

# Process the folder and create plots
process_folder(csv_folder_path, full_data, output_folder_path)

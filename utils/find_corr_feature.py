import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

full_data = pd.read_csv('data/svr/final_data.csv')

def process_folder(folder_path, full_data, output_folder, spearman_output_file):
    os.makedirs(output_folder, exist_ok=True)
    
    spearman_results = []

    for filename in os.listdir(folder_path):
        if filename.endswith('.csv'):
            print(filename)
            file_path = os.path.join(folder_path, filename)
            data = pd.read_csv(file_path)
            
            dockQ_values = []
            weighted_avg_values = []

            for _, row in data.iterrows():
                pdb_file = row['pdb_file']
                weighted_avg = row['patch_alignment_score']
                
                matching_row = full_data[full_data['pdb_file'] == pdb_file]
                
                if not matching_row.empty:
                    dockQ = matching_row['DockQ'].values[0]
                    dockQ_values.append(dockQ)
                    weighted_avg_values.append(weighted_avg)
    
            if dockQ_values and weighted_avg_values:
                spearman_corr, _ = spearmanr(weighted_avg_values, dockQ_values)
                spearman_results.append({'filename': filename, 'spearman_corr': spearman_corr})
            
            print("Plotting")
            plt.figure(figsize=(10, 6))
            plt.scatter(weighted_avg_values, dockQ_values, facecolors='none', edgecolors='black')
            plt.title(f'patch_alignment_score vs DockQ for {filename}')
            plt.xlabel('patch_alignment_score')
            plt.ylabel('DockQ')
            plt.grid(False)

            if weighted_avg_values:
                x_min, x_max = min(weighted_avg_values), max(weighted_avg_values)
                plt.xlim(x_min - 0.1*(x_max - x_min), x_max + 0.1*(x_max - x_min))
            if dockQ_values:
                y_min, y_max = min(dockQ_values), max(dockQ_values)
                plt.ylim(y_min - 0.1*(y_max - y_min), y_max + 0.1*(y_max - y_min))
            
            plot_filename = os.path.join(output_folder, f'{filename}.png')
            plt.savefig(plot_filename)
            print("Saved plot")
            plt.close()
            print()
    
    spearman_df = pd.DataFrame(spearman_results)
    spearman_df.to_csv(spearman_output_file, index=False)
    print(f"Spearman correlations saved to {spearman_output_file}")

csv_folder_path = 'all_csvs/patches'
output_folder_path = 'spearman_plots/patches/patch_alignment_score'
spearman_output_file = 'spearman_plots/patches/patch_alignment_score/patches_spearman.csv'

process_folder(csv_folder_path, full_data, output_folder_path, spearman_output_file)

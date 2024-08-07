import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

# Load the data from the CSV file
data = pd.read_csv('final_data_groups.csv')

# Group the data by 'pdb_id'
grouped = data.groupby('pdb_id')

# Initialize an empty list to store the results
correlation_results = []

# Create a scatter plot for each group
for pdb_id, group in grouped:
    rsm_values = group['rsm']
    dockq_values = group['DockQ']
    
    # Calculate the Spearman correlation
    correlation, p_value = spearmanr(rsm_values, dockq_values)
    correlation_results.append({'pdb_id': pdb_id, 'spearman_correlation': correlation})
    
    # Create scatter plot
    plt.figure(figsize=(8, 6))
    plt.scatter(rsm_values, dockq_values, edgecolors='black', facecolors='none')
    plt.title(f'{pdb_id} (Spearman Correlation: {correlation:.2f})', fontsize=20)
    plt.xlabel('RSM', fontsize=20)
    plt.ylabel('DockQ', fontsize=20)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

    plt.grid(False)
    
    # Save each plot as a file
    plt.savefig(f'/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/spearman_plots/rsm/{pdb_id}.png')
    plt.close()

# Convert the results to a DataFrame
results_df = pd.DataFrame(correlation_results)

# Display the results
print(results_df)

# Optionally, save the results to a CSV file
results_df.to_csv('rsm_spearman_correlations.csv', index=False)

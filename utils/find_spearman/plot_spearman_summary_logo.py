import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data from CSV files
methods = ['LOGO_C_R_F', 'LOGO_F', 'LOGO_R', 'LOGO_C']
data = {}

for method in methods:
    data[method] = pd.read_csv(f'all_predictions/{method}.csv')

# Combine data into a single DataFrame and calculate average Spearman correlation for each pdb_id
combined_df = pd.concat([data[method].set_index('pdb_id') for method in methods], axis=1, keys=methods)
combined_df.columns = combined_df.columns.map('_'.join)

# Multiply Spearman correlations by -1 if necessary
for method in methods:
    combined_df[f'{method}_spearman_correlation'] *= -1

# Load average Spearman values from another CSV
avg_spearman = pd.read_csv("naomi_spearman_dict.csv")

# Merge the average Spearman values into the combined DataFrame
for index, row in avg_spearman.iterrows():
    combined_df.at[row['pdb'], 'avg_spearman'] = row['avg_spearman']

# Handle missing values if any
combined_df['avg_spearman'] = combined_df['avg_spearman'].fillna(0)

# Sort by average Spearman correlation
df_sorted = combined_df.sort_values(by='avg_spearman', ascending=True).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(14, 7))

# Define markers and colors
markers = ['^', 'o', 's', 'P', '*', 'D', 'X', '<', '>']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple']

# Plot each method with different markers and colors
for i, method in enumerate(methods):
    ax.plot(df_sorted['pdb_id'], df_sorted[f'{method}_spearman_correlation'], 
            marker=markers[i % len(markers)], linestyle='None', 
            label=method, color=colors[i % len(colors)], 
            markersize=10, markerfacecolor='none', markeredgewidth=1.5)

# Plot the average Spearman correlation as a line
ax.plot(df_sorted['pdb_id'], df_sorted['avg_spearman'], linestyle='-', color='k', label='Average', linewidth=2)

# Customize the plot
ax.set_xlabel('CAPRI Targets')
ax.set_ylabel('Spearman Correlation')
ax.yaxis.grid(True)
ax.xaxis.grid(False)
ax.legend()
ax.set_ylim(-1, 0.4)
plt.xticks(rotation=45)

# Save the plot
plt.savefig('find_spearman/spearman_correlation_plot_logo.png')

# Show the plot
plt.show()

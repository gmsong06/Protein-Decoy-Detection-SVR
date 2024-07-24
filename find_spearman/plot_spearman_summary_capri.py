import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data from CSV files
# methods = ['CAPRI_C_R_F', '0.001_0.01', '0.001_0.001', '0.1_10', '0.001_10', '0.001_100']

methods = ['CAPRI_C_R_F']

gamma = [0.001, 0.01, 0.1, 1, 10, 100]
c = [0.001, 0.01, 0.1, 1, 10, 100]

for i in range(len(gamma)):
    for j in range(len(gamma)):
        methods.append(f"{gamma[i]}_{c[j]}")

data = {}
for method in methods:
    data[method] = pd.read_csv(f'all_predictions/{method}.csv')

# Combine data into a single DataFrame and calculate average Spearman correlation for each pdb_id
combined_df = pd.concat([data[method].set_index('pdb_id') for method in methods], axis=1, keys=methods)
combined_df.columns = combined_df.columns.map('_'.join)
for method in methods:
    combined_df[f'{method}_spearman_correlation'] *= -1  # Multiply Spearman correlations by -1
combined_df['average_spearman'] = combined_df[[f'{method}_spearman_correlation' for method in methods]].mean(axis=1)
combined_df['std_spearman'] = combined_df[[f'{method}_spearman_correlation' for method in methods]].std(axis=1)

# Sort by average Spearman correlation
df_sorted = combined_df.sort_values(by='average_spearman', ascending=True).reset_index()

# Create the plot
fig, ax = plt.subplots(figsize=(14, 7))

# Define markers and colors
markers = ['^', 'o', 's', 'P', '*', 'D', 'X', '<', '>', 'H', '1', '2', '3', '4']
colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'purple', 'brown', 'maroon', 'gold', 'silver', 'coral', 'turquoise', 'indigo', 'violet', 'crimson', 'orchid', 'sienna', 'plum', 'lavender']

# Plot each method with different markers and colors
for i, method in enumerate(methods):
    ax.plot(df_sorted['pdb_id'], df_sorted[f'{method}_spearman_correlation'], 
            marker=markers[i % len(markers)], linestyle='None', 
            label=method, color=colors[i % len(colors)], 
            markersize=10, markerfacecolor='none', markeredgewidth=1.5)

# Plot the average Spearman correlation with error bars
ax.errorbar(df_sorted['pdb_id'], df_sorted['average_spearman'], yerr=df_sorted['std_spearman'],
            fmt='', linestyle='-', color='k', label='Average', capsize=5)

# Customize the plot
ax.set_xlabel('CAPRI Targets')
ax.set_ylabel('Spearman')
ax.yaxis.grid(True)
ax.xaxis.grid(False)
ax.legend()
ax.set_ylim(-1, .4)

plt.xticks(rotation=45)

# Save the plot
plt.savefig('find_spearman/spearman_correlation_plot_logo.png')

# Show the plot
plt.show()
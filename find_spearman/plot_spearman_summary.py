import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data from CSV files
methods = ['CAPRI_ALL']
data = {}
for method in methods:
    data[method] = pd.read_csv(f'/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/all_predictions/{method}.csv')

# Extract pdb_ids and Spearman correlations
pdb_ids = data[methods[0]]['pdb_id']
spearman_correlations = data[methods[0]]['spearman_correlation']

# Calculate average Spearman correlation for each pdb_id (since only one method is provided, it's the same as spearman_correlations)
average_spearman = spearman_correlations

# Combine into a DataFrame
df = pd.DataFrame({'pdb_id': pdb_ids, 'spearman_correlation': spearman_correlations, 'average_spearman': average_spearman})

# Sort by average Spearman correlation
df_sorted = df.sort_values(by='average_spearman', ascending=False)

print(df_sorted)
# Create the plot
fig, ax = plt.subplots(figsize=(14, 7))

# Plot each method with different markers and colors
markers = ['^']
colors = ['b']

for i, method in enumerate(methods):
    ax.plot(df_sorted['pdb_id'], df_sorted['spearman_correlation'], marker=markers[i], linestyle='None', label=method, color=colors[i], markersize=10, markerfacecolor='none', markeredgewidth=1.5)

# Customize the plot
ax.set_ylabel('Spearman Correlation')
ax.yaxis.grid(True)
ax.xaxis.grid(False)
ax.legend()
ax.set_ylim(0, 1)

plt.xticks(rotation=45)

# Save the plot
plt.savefig('/home/annsong/Desktop/Yale_Research_Internship_24/Protein-Decoy-Detection-SVR/find_spearman/spearman_correlation_plot.png')

# Show the plot
plt.show()

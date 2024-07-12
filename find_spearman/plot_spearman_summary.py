import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

spearmans = pd.read_csv("spearman_plots/all_capri/spearman_correlations.csv")

methods = ['CAPRI_FULL']
keys = [f'Key {i+1}' for i in range(1)]


data = {method: np.random.uniform(-1, 1, 11) for method in methods}

for method in methods:
    data[method]

# Plot the data
fig, ax = plt.subplots(figsize=(14, 7))

markers = ['^', 'o', 's', 'P', '*', 'D']
colors = ['b', 'g', 'r', 'c', 'm', 'y']

for i, (method, spearman_values) in enumerate(data.items()):
    ax.plot(keys, spearman_values, marker=markers[i], linestyle='None', label=method, color=colors[i])

# Calculate and plot the average
average = np.mean(list(data.values()), axis=0)
ax.plot(keys, average, marker='x', linestyle='-', color='k', label='Average')

ax.set_xlabel('Keys')
ax.set_ylabel('Spearman Correlation (ρ)')
ax.set_title('Spearman Correlation for Different Methods')
ax.legend()
plt.xticks(rotation=45)
plt.grid(True)

plt.show()

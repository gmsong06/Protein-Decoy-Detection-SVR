# Normalizes data with z-score

import pandas as pd

df = pd.read_csv("data/final_data_capri.csv")

columns_to_normalize = ['contacts', 'rsm', 'flatness', 'hydro']

mean = df[columns_to_normalize].mean()
std = df[columns_to_normalize].std()

# Apply Z-score normalization to the specified columns
df[columns_to_normalize] = (df[columns_to_normalize] - mean) / std

# Display the updated DataFrame

df.to_csv("data/final_data_capri_normalized.csv", index=False)
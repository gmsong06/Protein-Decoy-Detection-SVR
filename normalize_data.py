import pandas as pd
import numpy as np

df = pd.read_csv("final_data_capri_with_hydro.csv")

columns_to_normalize = ['contacts', 'rsm', 'flatness', 'hydro']

mean = df[columns_to_normalize].mean()
std = df[columns_to_normalize].std()

# Apply Z-score normalization to the specified columns
df[columns_to_normalize] = (df[columns_to_normalize] - mean) / std

# Display the updated DataFrame

df.to_csv("final_data_capri_normalized.csv", index=False)
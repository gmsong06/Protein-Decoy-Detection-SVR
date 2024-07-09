import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--remove", nargs='+', type=str, help="Features to remove")
parser.add_argument("--output_path", type=str, help="Path to output prediction file")
args = parser.parse_args()

# Load the data
df = pd.read_csv('combined_data.csv')

specific_column = 'DockQ'
missing_in_specific_column = df[specific_column].isnull()

# Extracting the 'pdb_file' values where 'DockQ' is missing
missing_pdb_files = df[missing_in_specific_column]['pdb_file']

# Printing the missing 'pdb_file' values
print(missing_pdb_files)

# Saving the missing 'pdb_file' values to a text file
missing_pdb_files.to_csv('missing_pdb_files.txt', index=False, header=False)

to_remove = args.remove

if 'none' in to_remove:
    to_remove = []

to_remove.append('DockQ')
to_remove.append('pdb_file')
to_remove.append('pdb_id')

print(to_remove)
print(args.output_path)

X = df.drop(columns=to_remove)
y = df['DockQ']
groups = df['pdb_id']
pdb_files = df['pdb_file']

logo = LeaveOneGroupOut()

results = []
predictions = []

print("BEGINNING LOOPING")
for train_index, test_index in logo.split(X, y, groups):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    pdb_test = pdb_files.iloc[test_index]
    print("DATA SPLIT")
    
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('svr', SVR())
    ])
    
    print("BEGINNING TRAINING")
    pipeline.fit(X_train, y_train)
    print("DONE TRAINING")

    y_pred = pipeline.predict(X_test)
    print("DONE PREDICTING")
    
    mae = mean_absolute_error(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    results.append({
        'Group': groups.iloc[test_index].iloc[0],
        'MAE': mae,
        'MSE': mse,
        'R2': r2
    })
    
    for pdb, actual, pred in zip(pdb_test, y_test, y_pred):
        predictions.append({'pdb_file': pdb, 'actual_DockQ': actual, 'prediction': pred})

results_df = pd.DataFrame(results)
print(results_df)

predictions_df = pd.DataFrame(predictions)
predictions_df.to_csv(args.output_path, index=False)
print(predictions_df)

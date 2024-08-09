# Not in use

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.feature_selection import RFE
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import cross_val_score, learning_curve
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Command-line argument parsing
parser = argparse.ArgumentParser()
parser.add_argument("--remove", nargs='+', type=str, help="Features to remove")
parser.add_argument("--output_path", type=str, help="Path to output prediction file")
args = parser.parse_args()

# Load the data
train = pd.read_csv('final_data_groups_hydro.csv')
test = pd.read_csv('final_data_capri_with_hydro.csv')

specific_column = 'DockQ'
missing_in_specific_column = test[specific_column].isnull()

# Extracting the 'pdb_file' values where 'DockQ' is missing
missing_pdb_files = test[missing_in_specific_column]['pdb_file']
missing_pdb_files.to_csv('missing_pdb_files.txt', index=False, header=False)

to_remove = args.remove

if 'none' in to_remove:
    to_remove = []

to_remove.extend(['DockQ', 'pdb_file', 'pdb_id'])

X_train = train.drop(columns=to_remove)
y_train = train['DockQ']
X_test = test.drop(columns=to_remove)
y_test = test['DockQ']
pdb_files = test['pdb_file']

# Define feature ranking manually
ranking = [1, 1, 3, 2]  # Assuming you already have these rankings
weights = 1.0 / np.array(ranking)
weights = weights / np.sum(weights)

print(f"RANKING IS {ranking}")
print(f"WEIGHTS ARE {weights}")

# Apply weights to the features
X_train_weighted = X_train * weights
X_test_weighted = X_test * weights

# Model pipeline with SVR
svr = SVR(kernel='rbf', C=1.0)  # You can tune the C parameter
pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('svr', svr)
])

# Cross-validation scores
scores = cross_val_score(pipeline, X_train_weighted, y_train, cv=5, scoring='r2')
print(f"Cross-Validation R² scores: {scores}")
print(f"Mean CV R² score: {np.mean(scores)}")

# Train and evaluate the model
print("BEGINNING TRAINING")
pipeline.fit(X_train_weighted, y_train)
print("DONE TRAINING")

y_pred = pipeline.predict(X_test_weighted)
print("DONE PREDICTING")

mae = mean_absolute_error(y_test, y_pred)
mse = mean_squared_error(y_test, y_pred)
r2 = r2_score(y_test, y_pred)

print(f"Test MAE: {mae}")
print(f"Test MSE: {mse}")
print(f"Test R²: {r2}")

# Learning curves
train_sizes, train_scores, test_scores = learning_curve(pipeline, X_train_weighted, y_train, cv=5, scoring='r2', train_sizes=np.linspace(0.1, 1.0, 10))

train_scores_mean = np.mean(train_scores, axis=1)
test_scores_mean = np.mean(test_scores, axis=1)

plt.figure()
plt.plot(train_sizes, train_scores_mean, label='Training score')
plt.plot(train_sizes, test_scores_mean, label='Validation score')
plt.xlabel('Training Size')
plt.ylabel('R² Score')
plt.title('Learning Curves')
plt.legend()
plt.grid()
plt.show()

# Save predictions
predictions = pd.DataFrame({'pdb_file': pdb_files, 'actual_DockQ': y_test, 'prediction': y_pred})
predictions.to_csv(args.output_path, index=False)
print(predictions)

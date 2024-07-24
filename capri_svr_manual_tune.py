import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import argparse
import matplotlib.pyplot as plt

# Define the argument parser
parser = argparse.ArgumentParser()
parser.add_argument("--remove", nargs='+', type=str, help="Features to remove")
parser.add_argument("--output_path", type=str, help="Path to output prediction file")
parser.add_argument("--gamma", type=float, help="Gamma hyperparameter")
parser.add_argument("--c", type=float, help="C hyperparameter")
parser.add_argument("--fig_output_name", type=str, help="Name of figure output file")
args = parser.parse_args()

# Open a file to store print statements
log_file = open(f'{args.fig_output_name}.txt', 'w')

log_file.write("CREATED LOG FILE\n")

# Load the data
train = pd.read_csv('/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/final_data_groups_hydro.csv')
test = pd.read_csv('/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/final_data_capri_with_hydro.csv')

specific_column = 'DockQ'
missing_in_specific_column = test[specific_column].isnull()

# Extracting the 'pdb_file' values where 'DockQ' is missing
missing_pdb_files = test[missing_in_specific_column]['pdb_file']

# Logging the missing 'pdb_file' values
log_file.write("Missing pdb_file values:\n")
log_file.write(missing_pdb_files.to_string(index=False))
log_file.write("\n")
log_file.flush()

# Saving the missing 'pdb_file' values to a text file
missing_pdb_files.to_csv('missing_pdb_files.txt', index=False, header=False)

to_remove = args.remove

if 'none' in to_remove:
    to_remove = []

to_remove.append('DockQ')
to_remove.append('pdb_file')
to_remove.append('pdb_id')

log_file.write(f"Features to remove: {to_remove}\n")
log_file.write(f"Output path: {args.output_path}\n")
log_file.flush()

X_train = train.drop(columns=to_remove)
y_train = train['DockQ']
X_test = test.drop(columns=to_remove)
y_test = test['DockQ']
groups = train['pdb_id']
pdb_files = test['pdb_file']

results = []
predictions = []

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('svr', SVR(kernel='rbf', gamma=args.gamma, C=args.c))
])

# Train the model
pipeline.fit(X_train, y_train)

# Evaluate on training set
y_train_pred = pipeline.predict(X_train)
mae_train = mean_absolute_error(y_train, y_train_pred)
mse_train = mean_squared_error(y_train, y_train_pred)
r2_train = r2_score(y_train, y_train_pred)

# Evaluate on test set
y_pred = pipeline.predict(X_test)
log_file.write("DONE PREDICTING\n")
log_file.flush()  # Ensure immediate writing to the log file

mae_test = mean_absolute_error(y_test, y_pred)
mse_test = mean_squared_error(y_test, y_pred)
r2_test = r2_score(y_test, y_pred)

results.append({
    'MAE_Train': mae_train,
    'MSE_Train': mse_train,
    'R2_Train': r2_train,
    'MAE_Test': mae_test,
    'MSE_Test': mse_test,
    'R2_Test': r2_test,
    'Gamma': args.gamma,
    'C': args.c
})

for pdb, actual, pred in zip(pdb_files, y_test, y_pred):
    predictions.append({'pdb_file': pdb, 'actual_DockQ': actual, 'prediction': pred})

results_train = pd.DataFrame(results)
log_file.write(results_train.to_string(index=False))
log_file.write("\n")
log_file.flush()  # Ensure immediate writing to the log file

predictions_train = pd.DataFrame(predictions)
predictions_train.to_csv(args.output_path, index=False)
log_file.write(predictions_train.to_string(index=False))
log_file.write("\n")
log_file.flush()  # Ensure immediate writing to the log file

# Visualization
plt.figure(figsize=(14, 6))

# Training set results
plt.subplot(1, 2, 1)
plt.scatter(y_train, y_train_pred, edgecolors='black', facecolors='none', alpha=0.5)
plt.plot([y_train.min(), y_train.max()], [y_train.min(), y_train.max()], 'k--', lw=2)
plt.xlabel('Actual DockQ')
plt.ylabel('Predicted DockQ')
plt.title('Training Set')

# Test set results
plt.subplot(1, 2, 2)
plt.scatter(y_test, y_pred, edgecolors='black', facecolors='none', alpha=0.5)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2)
plt.xlabel('Actual DockQ')
plt.ylabel('Predicted DockQ')
plt.title('Test Set')

plt.tight_layout()
plt.savefig(f"{args.fig_output_name}.png")

# Close the log file
log_file.close()

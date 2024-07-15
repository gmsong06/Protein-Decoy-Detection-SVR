import pandas as pd
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
train = pd.read_csv('/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/final_data_groups.csv')
test = pd.read_csv('/home/as4643/palmer_scratch/Protein-Decoy-Detection-SVR/final_data_capri_groups.csv')

to_remove = args.remove

if 'none' in to_remove:
    to_remove = []

to_remove.append('DockQ')
to_remove.append('pdb_file')
to_remove.append('pdb_id')

print(to_remove)
print(args.output_path)

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
    ('svr', SVR(kernel='poly'))  # Example for polynomial kernel
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
    'MAE': mae,
    'MSE': mse,
    'R2': r2
})

for pdb, actual, pred in zip(pdb_files, y_test, y_pred):
    predictions.append({'pdb_file': pdb, 'actual_DockQ': actual, 'prediction': pred})

results_train = pd.DataFrame(results)
print(results_train)

predictions_train = pd.DataFrame(predictions)
predictions_train.to_csv(args.output_path, index=False)
print(predictions_train)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.svm import SVR
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# Load the data
df = pd.read_csv('final_data_groups.csv')

X = df.drop(columns=['DockQ', 'pdb_file', 'pdb_id'])
y = df['DockQ']
groups = df['pdb_id']

logo = LeaveOneGroupOut()

results = []

for train_index, test_index in logo.split(X, y, groups):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y.iloc[train_index], y.iloc[test_index]
    
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('svr', SVR())
    ])
    
    pipeline.fit(X_train, y_train)
    
    y_pred = pipeline.predict(X_test)
    
    mae = mean_absolute_error(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    
    results.append({
        'Group': groups.iloc[test_index].iloc[0],  # The group number for this test set
        'MAE': mae,
        'MSE': mse,
        'R2': r2
    })

# Convert results to a DataFrame for better readability
results_df = pd.DataFrame(results)
print(results_df)

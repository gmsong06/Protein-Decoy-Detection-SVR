import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# Load the data
df = pd.read_csv('final_combined_data.csv')

# Split the data
train, test = train_test_split(df, test_size=0.2, random_state=42)
X_train, X_test = train.drop(columns=['DockQ', 'pdb_file']).values, test.drop(columns=['DockQ', 'pdb_file']).values
y_train, y_test = train['DockQ'].values, test['DockQ'].values

missing_locations = df.isnull().stack()
missing_locations = missing_locations[missing_locations]
print(missing_locations)

# Standardize the data
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train the SVR model
svr_lin = SVR(kernel='linear')
svr_lin.fit(X_train_scaled, y_train)

# Predict on train and test sets
train['linear_svr_pred'] = svr_lin.predict(X_train_scaled)
test['linear_svr_pred'] = svr_lin.predict(X_test_scaled)

# Calculate residuals
train['residuals'] = train['DockQ'] - train['linear_svr_pred']
test['residuals'] = test['DockQ'] - test['linear_svr_pred']

# Visualize the performance

# Scatter plot of Predicted vs Actual values
plt.figure(figsize=(12, 6))
plt.scatter(train['DockQ'], train['linear_svr_pred'], label='Train', alpha=0.6)
plt.scatter(test['DockQ'], test['linear_svr_pred'], label='Test', alpha=0.6)
plt.plot([train['DockQ'].min(), train['DockQ'].max()], [train['DockQ'].min(), train['DockQ'].max()], 'k--', lw=2)
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.title('Predicted vs Actual Values')
plt.legend()
plt.show()

# Residual plot
plt.figure(figsize=(12, 6))
plt.scatter(train['linear_svr_pred'], train['residuals'], label='Train', alpha=0.6)
plt.scatter(test['linear_svr_pred'], test['residuals'], label='Test', alpha=0.6)
plt.axhline(y=0, color='r', linestyle='--')
plt.xlabel('Predicted Values')
plt.ylabel('Residuals')
plt.title('Residual Plot')
plt.legend()
plt.show()

# Histogram of residuals
plt.figure(figsize=(12, 6))
plt.hist(train['residuals'], bins=30, alpha=0.6, label='Train')
plt.hist(test['residuals'], bins=30, alpha=0.6, label='Test')
plt.xlabel('Residuals')
plt.ylabel('Frequency')
plt.title('Histogram of Residuals')
plt.legend()
plt.show()

# Calculate and print performance metrics
mae_train = mean_absolute_error(y_train, train['linear_svr_pred'])
mse_train = mean_squared_error(y_train, train['linear_svr_pred'])
r2_train = r2_score(y_train, train['linear_svr_pred'])

mae_test = mean_absolute_error(y_test, test['linear_svr_pred'])
mse_test = mean_squared_error(y_test, test['linear_svr_pred'])
r2_test = r2_score(y_test, test['linear_svr_pred'])

print(f'Train MAE: {mae_train}')
print(f'Train MSE: {mse_train}')
print(f'Train R²: {r2_train}')

print(f'Test MAE: {mae_test}')
print(f'Test MSE: {mse_test}')
print(f'Test R²: {r2_test}')

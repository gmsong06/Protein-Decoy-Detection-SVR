from sklearn.svm import SVR
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.decomposition import PCA
import numpy as np
import os
import pandas as pd

df = pd.read_csv('combined_data.csv')

train, test = train_test_split(df, test_size=0.2, random_state=42)

train = train.sort_values('all_contacts')
test = test.sort_values('all_contacts')

print(train.head())

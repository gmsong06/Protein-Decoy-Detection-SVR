from sklearn.svm import SVR
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
import numpy as np
import os
import pandas as pd

data = pd.read_csv('name.csv') #replace name with name of big csv file

X = data[['interface_flatness', 'surface_area', 'number_of_contacts']].values

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
pca = PCA(n_components=3)
principal_component = pca.fit_transform(X_scaled)

data['score'] = principal_component[:, 0]
data.to_csv('data_with_score.csv', index=False)

X = data[['interface_flatness', 'surface_area', 'number_of_contacts']].values
y = data['score'].ravel()


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

pipeline = Pipeline(steps=[('standardscaler', StandardScaler()),
                    ('svr', SVR(C= 1.0, epsilon=0.2))])

pipeline.fit(X_train, y_train)

y_pred = pipeline.predict(X_test)



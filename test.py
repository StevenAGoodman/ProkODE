import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt


# Define the model function
def model_func(x, *params):
    x = np.array(x).T
    params = np.array(params).T
    return x @ params

# Generate sample data
np.random.seed(0)  # For reproducibility
n_samples = 100000
n_features = 3

# Random independent variables
X = np.random.rand(n_samples, n_features) 

# True parameters
true_params = [2.,3.,4.]

# Generate dependent variable with some noise
y = []
for x in X:
    x = [float(n) for n in x] 
    y.append(model_func(x, 2.,3., 4.) + np.random.normal(0, 10))

p = [1] * n_features

# Fit the model
popt, pcov = sp.optimize.curve_fit(model_func, X.T, y, p0 = p)

# Display the results
print("Fitted parameters:", popt, "\ncov:", pcov)


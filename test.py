import numpy as np
import pandas as pd
import scipy as sp
import matplotlib.pyplot as plt


# Define the model function
def model_func(x, *params):
    return params[0] * x[0] + params[1] * x[1]

# Generate sample data
np.random.seed(0)  # For reproducibility
n_samples = 100
n_features = 2

# Random independent variables
X = np.random.rand(n_samples, n_features) * 10

# True parameters
true_params = np.array([])

# Generate dependent variable with some noise
y = model_func(X, 1.5, -2.0) + np.random.normal(0, .01, n_samples)

p = [1] * n_features

# Fit the model
popt, pcov = sp.optimize.curve_fit(model_func, X, y, p0 = p)

# Display the results
print("Fitted parameters:", popt)

# Optional: Plotting the results (only for 2D visualization if you choose 2 features)
if n_features == 2:
    plt.scatter(X[:, 1], y, label='Data')
    plt.plot(X[:, 1], model_func(X, *popt), color='red', label='Fitted Model')
    plt.xlabel('X1')
    plt.ylabel('y')
    plt.legend()
    plt.show()

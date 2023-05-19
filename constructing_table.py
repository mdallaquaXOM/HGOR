import numpy as np
from pyDOE2 import lhs

# range of variables
API = [38, 55]
gamma_g = [0.65, 1.]
temperature = [130, 300]

bounds = np.vstack((API, gamma_g, temperature)).T
X = lhs(3, samples=10, random_state=123)

for i in range(3):
    min_val = bounds[0, i]
    max_val = bounds[1, i]
    X[:, i] = X[:, i] * (max_val - min_val) + min_val



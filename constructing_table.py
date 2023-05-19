import numpy as np
from correlations.utils import sampling

# range of variables
API = [38, 55]
gamma_g = [0.65, 1.]
temperature = [130, 300]

bounds = np.vstack((API, gamma_g, temperature)).T

X = sampling(sampling_type='lhs', n=3, n_samples=10,
             random_state=123, iterations=None, bounds=bounds)

a = 0






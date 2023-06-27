import os
import pickle
import itertools
import numpy as np
from correlations.utils import excel_columns_map
from correlations.HGOR_script import PVTCORR_HGOR
from correlations.optimizer import optimizeParameter
from correlations.definitions import C_exp_rat_8_blasingame

# map of columns to read the Excel file
columnsNumbers, columnsNames = excel_columns_map()

source_opt = 'PVT_Data_updated'

pvtc = PVTCORR_HGOR(sat_pressure=None,
                    hgor=2500,
                    columns=columnsNumbers,
                    path='data',
                    files=[source_opt],
                    columnsNames=columnsNames,
                    categorical_columns=['well_name', 'formation', 'fluid'],
                    skiprows=0,
                    dataAugmentation=0)

# Optimizer
# 1 - Global optimization
# 2 - Trust-Region Constrained Algorithm (method='trust-constr')
# 3 - Sequential Least SQuares Programming (SLSQP) Algorithm (method='SLSQP')
# 4 - Unconstrained minimization: Broyden-Fletcher-Goldfarb-Shanno algorithm (method='BFGS')

print()
print('Exponential Rational 8 Optimization')
lb = [-1] * 8
ub = [1] * 8
bounds = [lb, ub]

C_new_8 = optimizeParameter(pvtc, opt_equation='Rog',
                            algorithm=1,
                            metric_func='LSE',
                            correlation_method={'principle': 'exponential_rational_8', 'variation': 'optimized'},
                            source=source_opt,
                            bounds=bounds,
                            correct_gamma=False,
                            verbose=True,
                            callback=False
                            )

new_parameters = {'Rog': {'exponential_rational_8': C_new_8.x,
                          }}

# is there a file with optimized values already?
path = r"optimizedParam/opt_results.pickle"
fileExist = os.path.isfile(path)

if fileExist:
    old_parameters = pickle.load(open(path, "rb"))
    old_parameters['Rog'] = new_parameters['Rog']
    new_parameters = old_parameters

pickle.dump(new_parameters, open(path, "wb"))
# np.save(r'optimizedParam/opt_results.npy',  new_parameter)

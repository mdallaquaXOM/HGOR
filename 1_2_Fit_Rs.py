import os
import pickle
import numpy as np

from correlations.HGOR_script import PVTCORR_HGOR
from correlations.optimizer import optimizeParameter
from correlations.definitions import C_exp_rat_8_blasingame

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rgo', 'psat', 'Bo', 'visc_o'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper'],
                    dataAugmentation=1)

source_opt = 'PVT_Data'
source_curve = 'PVT_Data'


# Optimizer
# 1 - Global optimization
# 2 - Trust-Region Constrained Algorithm (method='trust-constr')
# 3 - Sequential Least SQuares Programming (SLSQP) Algorithm (method='SLSQP')
# 4 - Unconstrained minimization: Broyden-Fletcher-Goldfarb-Shanno algorithm (method='BFGS')
print('Vasquez and Beggs Optimization')
C_new_VB = optimizeParameter(pvtc, opt_equation='Rgo',
                             algorithm=3,
                             metric_func='LSE',
                             correlation_method={'principle': 'vasquez_beggs', 'variation': 'optimized'},
                             source=source_opt,
                             bounds=([1e-2, 1., 15., 40.], [8e-2, 2., 30, 50]),
                             x_start=np.array([0.0178, 1.187, 23.931, 47.]))
print()
print('Exponential Rational 8 Optimization')
C_new_8 = optimizeParameter(pvtc, opt_equation='Rgo',
                            algorithm=4,
                            metric_func='LSE',
                            correlation_method={'principle': 'exponential_rational_8', 'variation': 'optimized'},
                            source=source_opt,
                            x_start=np.array(C_exp_rat_8_blasingame))

# print()
# print('Exponential Rational 16 Optimization')
# C_new_16 = optimizeParameter(pvtc, algorithm=1,
#                              metric_func='LSE',
#                              correlation_method={'principle': 'exponential_rational_16', 'variation': 'optimize'},
#                              source=source_opt,
#                             bounds=np.tile(np.array([[-4], [4]]), (1, 16)),
#                              x_start=np.array(C_exp_rat_16_blasingame))

new_parameters = {'Rgo': {'vasquez_beggs': C_new_VB.x,
                          'exponential_rational_8': C_new_8.x,
                          'exponential_rational_16': None,
                          'ace': None,
                          'datadriven': None,
                          }}


# is there a file with optimized values already?
path = r"optimizedParam/opt_results.pickle"
fileExist = os.path.isfile(path)


if fileExist:
    old_parameters = pickle.load(open(path, "rb"))
    old_parameters['Rgo'] = new_parameters['Rgo']
    new_parameters = old_parameters


pickle.dump(new_parameters, open(path, "wb"))
# np.save(r'optimizedParam/opt_results.npy',  new_parameter)

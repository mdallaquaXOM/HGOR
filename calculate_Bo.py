import pickle
from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter
from correlations.definitions import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat', 'Bo_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper'],
                    dataAugmentation=1)

source_opt = 'PVT_Data'
source_curve = 'PVT_Data'

# Perform EDA
# EDA_seaborn(pvtc.pvt_table)
# EDA_plotly(pvtc.pvt_table)


# Optimizer
# 1 - Global optimization
# 2 - Trust-Region Constrained Algorithm (method='trust-constr')
# 3 - Sequential Least SQuares Programming (SLSQP) Algorithm (method='SLSQP')
# 4 - Unconstrained minimization: Broyden-Fletcher-Goldfarb-Shanno algorithm (method='BFGS')
print('Vasquez and Beggs Optimization')
C_new_VB = optimizeParameter(pvtc, opt_equation='Rs',
                             algorithm=3,
                             metric_func='LSE',
                             correlation_method={'principle': 'vasquez_beggs', 'variation': 'optimize'},
                             source=source_opt,
                             bounds=([1e-2, 1., 15., 40.], [8e-2, 2., 30, 50]),
                             x_start=np.array([0.0178, 1.187, 23.931, 47.]))
print()
print('Exponential Rational 8 Optimization')
C_new_8 = optimizeParameter(pvtc, opt_equation='Rs',
                            algorithm=4,
                            metric_func='LSE',
                            correlation_method={'principle': 'exponential_rational_8', 'variation': 'optimize'},
                            source=source_opt,
                            x_start=np.array(C_exp_rat_8_blasingame))

# print()
# print('Exponential Rational 16 Optimization')
# C_new_16 = optimizeParameter(pvtc, algorithm=1,
#                              metric_func='LSE',
#                              correlation_method={'principle': 'exponential_rational_16', 'variation': 'optimize'},
#                              source=source_opt,
#                              bounds=np.tile(np.array([[-4], [4]]), (1, 16)),
#                              x_start=np.array(C_exp_rat_16_blasingame))

new_parameter = {'Vasquez_Beggs': C_new_VB.x,
                 'exponential_rational_8': C_new_8.x,
                 # 'exponential_rational_16': C_new_16.x
                 }

# Saving values

pickle.dump(new_parameter, open(r"optimizedParam/opt_results.pickle", "wb"))
# np.save(r'optimizedParam/opt_results.npy',  new_parameter)



# Calculate RS
Rs, Rs_metrics = pvtc.compute_RS_values(new_parameter, source=source_curve)

# plots
plot_log_log(Rs, measured='Rs',
             calculated=['VB_original',
                         'VB_optimized',
                         'Exp_Rational_8_paper',
                         'Exp_Rational_8_optimized',
                         'Exp_Rational_16_paper',
                         'Exp_Rational_16_ed'
                         # 'Exp_Rational_16_optimized'
                         ],
             metrics_df=Rs_metrics,
             title='Rs (scf/stb) at saturation pressure')

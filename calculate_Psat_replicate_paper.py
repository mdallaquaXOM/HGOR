from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter
from correlations.definitions import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat', 'Bo_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper'],
                    dataAugmentation=0)

# optimize parameters for exponential Rational 16
# print()
# print('Exponential Rational 16 Optimization...',end="")
# C_new_16 = optimizeParameter(pvtc, algorithm=2,
#                              opt_equation='pb',
#                              metric_func='LSE',
#                              correlation_method='exponential_rational_16_optimize',
#                              source='PVT_paper',
#                              bounds=np.tile(np.array([[-3], [3]]), (1, 16)),
#                              x_start=C_exp_rat_16_blasingame)
print('Done!')
# new_parameter = {
#     # 'Vasquez_Beggs': C_new_VB.x,
#     'exponential_rational_16': C_new_16.x,
# }

# Calculate bubble point
pb, pb_metrics = pvtc.compute_pb_values(source='PVT_Data')

plot_properties(pb,
                measured='pressure',
                calculated=[
                    'VB_original',
                    'VB_paper',
                    'Exp_Rational_8_paper',
                    'Exp_Rational_16_paper',
                    'Exp_Rational_16_Ed'
                ],
                metrics_df=pb_metrics,
                title='Saturation Point [psia]',
                property='p_sat')

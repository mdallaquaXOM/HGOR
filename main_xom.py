from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat', 'Bo_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper'])

# Perform EDA
# EDA_seaborn(pvtc.pvt_table)
# EDA_plotly(pvtc.pvt_table)


# Optimizer
# 1 - Trust-Region Constrained Algorithm (method='trust-constr')
# 2 - Sequential Least SQuares Programming (SLSQP) Algorithm (method='SLSQP')
# 3 - Global optimization
C_new_vasquez = optimizeParameter(pvtc, type=1, metric_func='LSE', source='PVT_paper')

new_parameter = {'Vasquez_Beggs': C_new_vasquez.x}

# Calculate RS
Rs, Rs_metrics = pvtc.compute_RS_values(new_parameter, source='PVT_paper')


plot_log_log(Rs, measured='Rs',
             calculated=['VB_original',
                         'VB_modified',
                         'VB_paper',
                         'Exp_Rational_8',
                         'Exp_Rational_16_paper',
                         'Exp_Rational_16_new'],
             metrics_df=Rs_metrics,
             title='Rs (scf/stb) at saturation pressure')

plot_log_log(Rs, measured='Rs',
             calculated=['VB_original',
                         'VB_modified',
                         'Exp_Rational_8'],
             metrics_df=Rs_metrics,
             title='Rs (scf/stb) at saturation pressure')



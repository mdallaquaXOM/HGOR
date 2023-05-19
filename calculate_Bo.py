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


correlations = {'Bo': [{'principle': 'vasquez_beggs', 'variation': 'original'},
                       {'principle': 'vasquez_beggs', 'variation': 'rs_update'},
                       {'principle': 'exponential_rational_15', 'variation': 'michael'}],
                }

# Calculate RS
Rs, Rs_metrics = pvtc.compute_Bo_values(correlations, source=source_curve)

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

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


properties = {'Bo': [{'principle': 'vasquez_beggs', 'variation': 'original'},
                     {'principle': 'vasquez_beggs', 'variation': 'rs_update'},
                     {'principle': 'exponential_rational_15', 'variation': 'michael'}],
              }

# Calculate Bo
pvt_prop, pvt_metrics = pvtc.compute_PVT_Correlations(properties,
                                                      source=source_curve,
                                                      rs_best_correlation={'principle': 'exponential_rational_8',
                                                                           'variation': 'optimized'})

colums2plot = pvt_prop['Bo'].drop(['measured', 'HGOR'], axis=1).columns.values

plot_log_log(pvt_prop['Bo'], measured='measured',
             calculated=colums2plot,
             metrics_df=pvt_metrics['Bo'],
             title='Bo (Rb/stb) at saturation pressure',
             log_axis=False)

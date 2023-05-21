import pickle
from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter
from correlations.definitions import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat', 'Bo_psat', 'Bg_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper']
                    )

source_plots = 'PVT_Data'

# Correlations
correlations = {'Rs': [{'principle': 'vasquez_beggs', 'variation': 'original'},
                       {'principle': 'vasquez_beggs', 'variation': 'optimized'},
                       {'principle': 'vasquez_beggs', 'variation': 'paper'},
                       {'principle': 'exponential_rational_8', 'variation': 'meija'},
                       {'principle': 'exponential_rational_8', 'variation': 'optimized'},
                       {'principle': 'exponential_rational_16', 'variation': 'blasingame'},
                       {'principle': 'exponential_rational_16', 'variation': 'michael'}],
                }



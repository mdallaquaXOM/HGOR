from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter
from correlations.definitions import *

# from page 108


# range of variables
p_bounds = [4000., 7500]

bounds = {'temperature': [75, 225],
          'gamma_s': [0.65, 0.85],
          'API': [40, 60],
          }
# 'Rvi': [6.66667E-05, 0.00025]

# samples of API, Specific Gravity and Temp
inputs, _ = sampling(sampling_type='lhs', n_samples=len(bounds) * 10,
                     psat_bounds=p_bounds, n_psat=15,
                     random_state=123, bounds=bounds)

# New Correlations
properties = {
    # 'Rs': [
    #     {'principle': 'nasser', 'variation': 'GC'},
    #     {'principle': 'nasser', 'variation': 'VO'},
    # ],
    # 'Bo': [
    #     {'principle': 'nasser', 'variation': 'GC'},
    #     {'principle': 'nasser', 'variation': 'VO'},
    # ],
    'CGR': [
        # {'principle': 'nasser', 'variation': 'knownPsat'},
        {'principle': 'ace', 'variation': 'ovalle'},
    ],
    # 'Bg': [
    #     {'principle': 'nasser', 'variation': 'GC'},
    #     {'principle': 'nasser', 'variation': 'VO'},
    # ],
}

# pvtc = PVTCORR_HGOR(sat_pressure=None,
#                     Tsp=110, Psp=800,
#                     Tsp2=75, Psp2=60,
#                     hgor=2000,
#                     path='data',
#                     files=['PVT_Book'],
#                     correct_specific_gravity=False)

pvtc = PVTCORR_HGOR(sat_pressure=None,
                    Tsp=110, Psp=800,
                    hgor=2000,
                    data=inputs,
                    correct_specific_gravity=False)
# NEW Properties
pvtc.pvt_table['gamma_c'] = pvtc.pvt_table['gamma_s']
pvt_df = pvtc.compute_PVT_Correlations(properties)

input_table = pvtc.pvt_table

for property_, correlations in pvt_df.items():
    plot_synthetic_data(correlations, input_table, name=property_)

    plot_synthetic_data(correlations, input_table, name=property_, hueplot='sample')

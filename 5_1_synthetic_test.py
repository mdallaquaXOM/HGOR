import numpy as np
from correlations.utils import sampling_old, sampling
from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import plot_synthetic_data

# range of variables
# p_bounds = [100., 3500.]
#
# bounds = {'API': [38, 55.],
#           'gamma_s': [0.65, 1.],
#           'temperature': [130, 300],
#           }

p_bounds = [100., 9000.]

bounds = {'API': [38, 55.],
          'gamma_s': [0.65, 1.],
          'temperature': [130, 300],
          }


# samples of API, Specific Gravity and Temp
# inputs, _ = sampling_old(sampling_type='lhs', nVariables=3, n_samples=10, n_psat=100,
#                          random_state=123, bounds=bounds)

# samples of API, Specific Gravity and Temp
inputs, _ = sampling(sampling_type='lhs', n_samples=10,
                     psat_bounds=p_bounds, n_psat=15,
                     random_state=123, bounds=bounds)

# create class
pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    data=inputs
                    )

properties = {
    'Rgo': [
            {'principle': 'ace', 'variation': 'None'},
            {'principle': 'datadriven', 'variation': 'ann'},
            {'principle': 'datadriven', 'variation': 'randomforest'},
            {'principle': 'exponential_rational_8', 'variation': 'optimized'},
            # {'principle': 'exponential_rational_8', 'variation': 'blasingame'},
            # {'principle': 'exponential_rational_16', 'variation': 'blasingame'},
            # {'principle': 'exponential_rational_16', 'variation': 'michael'},
        ],
    'Bo': [
            # {'principle': 'vasquez_beggs', 'variation': 'original'},
            {'principle': 'vasquez_beggs', 'variation': 'rs_update'},
            {'principle': 'exponential_rational_15', 'variation': 'michael'}
        ],
    'visc_o': [
            {'principle': 'Beggs_and_Robinson', 'variation': 'original'},
            {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
            {'principle': 'exponential_rational_15', 'variation': 'michael'}
        ],
    'Rog': [
        {'principle': 'exponential_rational_8', 'variation': 'optimized'},
        ]
    }

pvt_df = pvtc.compute_PVT_Correlations(properties, rs_best_correlation={'principle': 'exponential_rational_8',
                                                                        'variation': 'optimized'})
#
# for property_, correlations in pvt_df.items():
#     plot_synthetic_data(correlations, inputs, name=property_,
#                         jumpLog='exponential_rational_16_michael', hueplot='sample')

for property_, correlations in pvt_df.items():
    print(property_)

    plot_synthetic_data(correlations, inputs, name=property_,
                        jumpLog='exponential_rational_16_michael')
    plot_synthetic_data(correlations, inputs, name=property_, jumpLog='exponential_rational_16_michael',
                        hueplot='sample')

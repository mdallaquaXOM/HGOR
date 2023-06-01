import numpy as np
from correlations.utils import sampling
from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import plot_synthetic_data

# range of variables
bounds = {'p': [100., 3500.],
          'API': [38, 55.],
          'gamma_s': [0.65, 1.],
          'temperature': [130, 300],
          }

# samples of API, Specific Gravity and Temp
inputs = sampling(sampling_type='lhs', nVariables=3, n_samples=10, n_psat=100,
                  random_state=123, bounds=bounds)

# create class
pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                     hgor=2000,
                     data=inputs
                     )

properties = {'Rs': [
    {'principle': 'vasquez_beggs', 'variation': 'original'},
    {'principle': 'exponential_rational_8', 'variation': 'blasingame'},
    {'principle': 'exponential_rational_8', 'variation': 'optimized'},
    {'principle': 'exponential_rational_16', 'variation': 'blasingame'},
    {'principle': 'exponential_rational_16', 'variation': 'michael'},
],
    'Bo': [
        {'principle': 'vasquez_beggs', 'variation': 'original'},
        {'principle': 'vasquez_beggs', 'variation': 'rs_update'},
        {'principle': 'exponential_rational_15', 'variation': 'michael'}
    ],
    'muob': [
        {'principle': 'Beggs_and_Robinson', 'variation': 'original'},
        {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
        {'principle': 'exponential_rational_15', 'variation': 'michael'}],
}

pvt_df = pvtc.compute_PVT_Correlations_v2(properties, rs_best_correlation={'principle': 'exponential_rational_8',
                                                                         'variation': 'optimized'})
#
# for property_, correlations in pvt_df.items():
#     plot_synthetic_data(correlations, inputs, name=property_,
#                         jumpLog='exponential_rational_16_michael', hueplot='sample')

for property_, correlations in pvt_df.items():
    plot_synthetic_data(correlations, inputs, name=property_,
                        jumpLog='exponential_rational_16_michael')

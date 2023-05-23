from correlations.utils import sampling, printInputValues, plot_comparePVT
from correlations.new_correlations import PVTCORR_HGOR

# New Correlations  - select JUST ONE !!!!!!
properties = {
    'Rs':
    # {'principle': 'vasquez_beggs', 'variation': 'original'},
    # {'principle': 'exponential_rational_8', 'variation': 'blasingame'},
        {'principle': 'exponential_rational_8', 'variation': 'optimized'},
    # {'principle': 'exponential_rational_16', 'variation': 'blasingame'},
    # {'principle': 'exponential_rational_16', 'variation': 'michael'},
    'Bo':
    # {'principle': 'vasquez_beggs', 'variation': 'original'},
        {'principle': 'vasquez_beggs', 'variation': 'rs_update'},
    # {'principle': 'exponential_rational_15', 'variation': 'michael'}
    'muob':
    # {'principle': 'Beggs_and_Robinson', 'variation': 'original'},
    {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
    #     {'principle': 'exponential_rational_15', 'variation': 'michael'}
}

# New Correlation class

bounds = {'p_sat': [100., 3500.],
          'API': [38, 55.],
          'gamma_s': [0.65, 1.],
          'temperature': [130, 300],
          }

# samples of API, Specific Gravity and Temp
inputs = sampling(sampling_type='lhs', nVariables=3, n_samples=1, n_psat=50,
                  random_state=538, bounds=bounds)

# create class
pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000.,
                    data=inputs
                    )

# Old Properties
pvt_old = pvtc.construct_PVT_table_old()

# Optimizer definitions
ranges = [(30, 55), (0.65, 1.2), (130, 300)]

input_star = pvtc.match_PVT_valuesHGOR(ranges,
                                       pvt_old=pvt_old,
                                       properties=properties,
                                       additional_details=True)

# Print comparison dor (API, Specific_Gravity, Temperature)
printInputValues(old=inputs, new=input_star)

# New Properties
pvt_new = pvtc.construct_PVT_table_new(properties, inputValues=input_star)

# Comparing both plots.
plot_comparePVT(inputs=pvtc.pvt_table,
                df_old=pvt_old,
                df_new=pvt_new,
                title='Match')


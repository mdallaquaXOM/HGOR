import pandas as pd
from correlations.utils import sampling
from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT

# range of variables
bounds = {'p_sat': [100., 3500.],
          'API': [38, 55.],
          'gamma_s': [0.65, 1.],
          'temperature': [130, 300],
          }

# samples of API, Specific Gravity and Temp
inputs = sampling(sampling_type='lhs', nVariables=3, n_samples=1, n_psat=100,
                  random_state=123, bounds=bounds)

# create class
pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                     hgor=2000,
                     data=inputs
                     )


# pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
#                     hgor=2000,
#                     columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat',
#                              'Bo_psat', 'Bg_psat', 'visc_o_psat'],
#                     path='data',
#                     files=['PVT_Data'],
#                     dataAugmentation=0)


# Old Properties
pvt_old = pvtc.construct_PVT_table_old()

# New Correlations  - JUST ONE !!!!!!
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
    # {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
        {'principle': 'exponential_rational_15', 'variation': 'michael'}
}

pvt_new = pvtc.construct_PVT_table_new(properties)

# Calculate metrics


# Comparing both plots.
plot_comparePVT(inputs=pvtc.pvt_table,
                df_old=pvt_old,
                df_new=pvt_new)

# include pressure
psat = pvtc.pvt_table['p_sat']
pvt_old.insert(0, 'p_sat', psat)
pvt_new.insert(0, 'p_sat', psat)

pvt_new.to_csv('outputs/pvt_new.csv', index=False)
pvt_old.to_csv('outputs/pvt_old.csv', index=False)
# pvtc.pvt_table.to_csv('outputs/inputs.csv', index=False)
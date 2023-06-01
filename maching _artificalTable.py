import pandas as pd

from correlations.utils import sampling, printInputValues, plot_comparePVT, metrics, metric2df
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

bounds = {'p': [100., 3500.],
          'API': [38, 55.],
          'gamma_s': [0.65, 1.],
          'temperature': [130, 300],
          }

# rukle of thumb: a sample size of 10 times the number of design variables
# https://designbuilder.co.uk/helpv7.0/Content/UASACalculationOptionsGeneral.htm#:~:text=4%2DLHS%3A%20Latin
# %20hypercube%20sampling,value%20of%20requested%20distribution%20range.
n_samples = 30

# samples of API, Specific Gravity and Temp
inputs, sampled_valeus = sampling(sampling_type='lhs', nVariables=3, n_samples=n_samples, n_psat=25,
                                  random_state=538, bounds=bounds)

# Optimizer definitions
ranges = [(30, 55), (0.65, 1.2), (130, 300)]

errors = []

for n_sample in range(n_samples):
    # create class
    pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                        hgor=2000.,
                        data=inputs[inputs['sample'] == n_sample]
                        )

    # Old Properties
    pvt_old = pvtc.construct_PVT_table_old()

    columnToMatch = ['Rgo', 'Bo']
    input_star, error = pvtc.match_PVT_valuesHGOR(ranges,
                                                  pvt_old=pvt_old,
                                                  properties=properties,
                                                  additional_details=True,
                                                  columnToMatch=columnToMatch
                                                  )

    # Print comparison dor (API, Specific_Gravity, Temperature)
    printInputValues(old=inputs, new=input_star)

    # New Properties
    pvt_new = pvtc.construct_PVT_table_new(properties, inputValues=input_star)

    # Get the error
    metrics_dict = metrics(pvt_old, pvt_new)
    metric_df = metric2df(metrics_dict, error)

    # Comparing both plots.
    plot_comparePVT(inputs=pvtc.pvt_table,
                    df_old=pvt_old,
                    df_new=pvt_new,
                    title=f'Match_sample_{n_sample}')

    # save errors
    errors.append(metric_df)

# What were the sample values
metric_all = pd.concat(errors).reset_index(drop=True)
summary = pd.concat([sampled_valeus, metric_all], axis=1)

summary.to_excel(fr"outputs/errors_match.xlsx")
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

bounds = {'p': [100., 4900.],
          'API': [38, 55.],
          'gamma_s': [0.65, 1.],
          'temperature': [130, 300],
          }

# rule of thumb: a sample size of 10 times the number of design variables
# https://designbuilder.co.uk/helpv7.0/Content/UASACalculationOptionsGeneral.htm#:~:text=4%2DLHS%3A%20Latin
# %20hypercube%20sampling,value%20of%20requested%20distribution%20range.
n_samples = 30

# samples of API, Specific Gravity and Temp
inputs, sampled_values = sampling(sampling_type='lhs', nVariables=3, n_samples=n_samples, n_psat=25,
                                  random_state=123, bounds=bounds)

# Optimizer definitions
# ranges = [(30, 55), (0.65, 1.2), (130, 300)]
ranges = [(20, 55), (0.5, 1.5), (100, 350)]

errors = []
matched_values = []

for n_sample in range(n_samples):
    # create class
    pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                        hgor=2000.,
                        data=inputs[inputs['sample'] == n_sample]
                        )
    pvt_input = pvtc.pvt_table

    # Old Properties
    pvt_old = pvtc.construct_PVT_table_old()

    # NEW Properties
    pvt_new = pvtc.construct_PVT_table_new(properties)

    # Matched PVT
    # print(sampled_values.iloc[n_sample, :].to_numpy())

    columnToMatch = ['Rgo', 'Bo']
    input_star, error = pvtc.match_PVT_valuesHGOR(ranges,
                                                  pvt_old=pvt_old,
                                                  properties=properties,
                                                  additional_details=True,
                                                  columnToMatch=columnToMatch,
                                                  disp=True,
                                                  x_start=sampled_values.iloc[n_sample, :].to_numpy(),
                                                  printXk=False
                                                  )
    matched_values.append(input_star)

    # Print comparison dor (API, Specific_Gravity, Temperature)
    printInputValues(old=inputs, new=input_star)

    # New Properties
    pvt_match = pvtc.construct_PVT_table_new(properties, inputValues=input_star)

    # Get the error
    metrics_dict = metrics(pvt_old, pvt_match)
    metric_df = metric2df(metrics_dict, error)

    # Comparing both plots.
    plot_comparePVT(inputs=pvt_input,
                    df_old=pvt_old,
                    df_new=pvt_new,
                    df_opt=pvt_match,
                    title=f'Match_sample_{n_sample}',
                    path=r'figures/matching/')

    # save errors
    errors.append(metric_df)

    # print tables
    psat = pvtc.pvt_table['p']
    pvt_old.insert(0, 'p', psat)
    pvt_match.insert(0, 'p', psat)

    with pd.ExcelWriter(fr'outputs/samples/pvt_{n_sample}.xlsx') as writer:
        pvt_old.to_excel(writer, index=False, sheet_name='original')
        pvt_match.to_excel(writer, index=False, sheet_name='matched')

# What were the sampled values
matched_values_df = pd.DataFrame(matched_values, columns=sampled_values.columns)

# add second layer of index
sampled_values.columns = pd.MultiIndex.from_product([['original'], sampled_values.columns])
matched_values_df.columns = pd.MultiIndex.from_product([['matched'], matched_values_df.columns])

metric_all = pd.concat(errors).reset_index(drop=True)
summary = pd.concat([sampled_values, matched_values_df, metric_all], axis=1)

summary.to_excel(fr"outputs/errors_match.xlsx")

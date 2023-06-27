import pandas as pd

from correlations.utils import sampling_old, printInputValues, plot_comparePVT, metrics, metric2df, relativeErrorforMatch
from correlations.HGOR_script import PVTCORR_HGOR

# New Correlations  - select JUST ONE !!!!!!
properties = {
    'Rgo':
    # {'principle': 'vasquez_beggs', 'variation': 'original'},
    # {'principle': 'exponential_rational_8', 'variation': 'blasingame'},
        {'principle': 'exponential_rational_8', 'variation': 'optimized'},
    # {'principle': 'exponential_rational_16', 'variation': 'blasingame'},
    # {'principle': 'exponential_rational_16', 'variation': 'michael'},
    'Bo':
    # {'principle': 'vasquez_beggs', 'variation': 'original'},
        {'principle': 'vasquez_beggs', 'variation': 'rs_update'},
    # {'principle': 'exponential_rational_15', 'variation': 'michael'}
    'visc_o':
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
inputs, input_sampled = sampling_old(sampling_type='lhs', nVariables=3, n_samples=n_samples, n_psat=15,
                                     random_state=123, bounds=bounds)

# Optimizer definitions
# ranges = [(30, 55), (0.65, 1.2), (130, 300)]
ranges = [(20, 55), (0.45, 1.8), (100, 400)]

errors = []
input_matched = []
columnToMatch = ['Rgo', 'Bg']

for n_sample in range(n_samples):
    # create class
    pvt_input = inputs[inputs['sample'] == n_sample].reset_index(drop=True)

    pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                        hgor=2000.,
                        data=pvt_input
                        )

    # Old Properties
    pvt_old = pvtc.construct_PVT_table_old()

    # NEW Properties
    pvt_new = pvtc.construct_PVT_table_new(properties)

    # Matched PVT
    # print(sampled_values.iloc[n_sample, :].to_numpy())

    input_star, error_opt = pvtc.match_PVT_valuesHGOR(ranges,
                                                      pvt_old=pvt_old,
                                                      properties=properties,
                                                      additional_details=True,
                                                      columnToMatch=columnToMatch,
                                                      disp=True,
                                                      x_start=input_sampled.iloc[n_sample, :].to_numpy(),
                                                      printXk=False
                                                      )
    input_matched.append(input_star)

    # Print comparison dor (API, Specific_Gravity, Temperature)
    printInputValues(old=inputs, new=input_star)

    # New Properties
    pvt_match = pvtc.construct_PVT_table_new(properties, inputValues=input_star)

    # Get the error and metrics
    metrics_dict = metrics(pvt_old, pvt_match)
    metric_df = metric2df(metrics_dict).reset_index(drop=True)

    error_total, _ = relativeErrorforMatch(pvt_old, pvt_match)
    error_df = pd.DataFrame([[error_total, error_opt]], columns=['total', 'optimization'])
    error_df.columns = pd.MultiIndex.from_product([['errors'], error_df.columns])

    errors.append(pd.concat([error_df, metric_df], axis=1))

    # Comparing both plots.
    plot_comparePVT(inputs=pvt_input,
                    df_old=pvt_old,
                    df_new=pvt_new,
                    df_opt=pvt_match,
                    title=f'Match_sample_{n_sample}',
                    path=r'figures/matching/')

    # print tables
    psat = pvtc.pvt_table['p']
    pvt_old.insert(0, 'psat', psat)
    pvt_match.insert(0, 'psat', psat)

    with pd.ExcelWriter(fr'outputs/samples/pvt_{n_sample}.xlsx') as writer:
        pvt_old.to_excel(writer, index=False, sheet_name='original')
        pvt_match.to_excel(writer, index=False, sheet_name='matched')

# What were the sampled values
matched_values_df = pd.DataFrame(input_matched, columns=input_sampled.columns)

# add second layer of index
input_sampled.columns = pd.MultiIndex.from_product([['original'], input_sampled.columns])
matched_values_df.columns = pd.MultiIndex.from_product([['matched'], matched_values_df.columns])

metric_all = pd.concat(errors).reset_index(drop=True)
summary = pd.concat([input_sampled, matched_values_df, metric_all], axis=1)

summary.to_excel(fr"outputs/errors_match.xlsx")

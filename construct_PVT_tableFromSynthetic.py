from correlations.utils import sampling
from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT, metrics, metric2df
import pandas as pd

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
        # {'principle': 'exponential_rational_15', 'variation': 'michael'}
}

# range of variables
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
inputs, sampled_values = sampling(sampling_type='lhs', nVariables=3, n_samples=n_samples, n_psat=50,
                                  random_state=123, bounds=bounds)

errors = []
for n_sample in range(n_samples):
    # overwrite inputs
    # inputs['API'] = 51.76405
    # inputs['temperature'] = 261.4015
    # inputs['gamma_s'] = 0.7932375

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

    # Comparing both plots.
    plot_comparePVT(inputs=pvt_input,
                    df_old=pvt_old,
                    df_new=pvt_new,
                    path=r'figures/synthetic/',
                    title=f'Synthetic {n_sample}')

    # Get the error
    metrics_dict = metrics(pvt_old, pvt_new)
    metric_df = metric2df(metrics_dict)

    errors.append(metric_df)

    # include pressure
    # psat = pvtc.pvt_table['p_sat']
    # pvt_old.insert(0, 'p_sat', psat)
    # pvt_new.insert(0, 'p_sat', psat)
    #
    # pvt_new.to_csv(fr'outputs/pvt_new_{type_of_data}.csv', index=False)
    # pvt_old.to_csv(fr'outputs/pvt_old{type_of_data}.csv', index=False)
    # pvtc.pvt_table.to_csv(fr'outputs/inputs_{type_of_data}.csv', index=False)

    del pvtc

metric_all = pd.concat(errors).reset_index(drop=True)
summary = pd.concat([sampled_values, metric_all], axis=1)

summary.to_excel(fr"outputs/errors_old_new.xlsx")

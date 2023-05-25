from correlations.utils import sampling
from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT

type_of_datas = ['Synthetic', 'Lab']

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
    # {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
        {'principle': 'exponential_rational_15', 'variation': 'michael'}
}

Tsp = 14.7
Psp = 60
hgor_treshold = 2000

for i_type, type_of_data in enumerate(type_of_datas):

    if type_of_data.lower() == 'synthetic':
        # range of variables
        bounds = {'p': [100., 3500.],
                  'API': [38, 55.],
                  'gamma_s': [0.65, 1.],
                  'temperature': [130, 300],
                  }

        # samples of API, Specific Gravity and Temp
        inputs, _ = sampling(sampling_type='lhs', nVariables=3, n_samples=1, n_psat=50,
                             random_state=123, bounds=bounds)

        # overwrite inputs
        inputs['API'] = 51.76405
        inputs['temperature'] = 261.4015
        inputs['gamma_s'] = 0.7932375

        # create class
        pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=Tsp, Psp=Psp,
                            hgor=hgor_treshold,
                            data=inputs
                            )
    else:
        pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=Tsp, Psp=Psp,
                            hgor=hgor_treshold,
                            columns=['temperature', 'API', 'gamma_s', 'Rgo', 'p',
                                     'Bo', 'Bg', 'visc_o'],
                            path='data',
                            files=['PVT_Data'],
                            dataAugmentation=0)

    # Old Properties
    pvt_old = pvtc.construct_PVT_table_old()

    # NEW Properties
    pvt_new = pvtc.construct_PVT_table_new(properties)

    # Calculate metrics

    # Comparing both plots.
    plot_comparePVT(inputs=pvtc.pvt_table,
                    df_old=pvt_old,
                    df_new=pvt_new,
                    title=type_of_data)

    # include pressure
    # psat = pvtc.pvt_table['p_sat']
    # pvt_old.insert(0, 'p_sat', psat)
    # pvt_new.insert(0, 'p_sat', psat)
    #
    # pvt_new.to_csv(fr'outputs/pvt_new_{type_of_data}.csv', index=False)
    # pvt_old.to_csv(fr'outputs/pvt_old{type_of_data}.csv', index=False)
    # pvtc.pvt_table.to_csv(fr'outputs/inputs_{type_of_data}.csv', index=False)

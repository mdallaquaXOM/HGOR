from correlations.utils import sampling
from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT

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
plot_comparePVT(properties=['Rgo', 'Bo', 'visc_o'],
                inputs=pvtc.pvt_table,
                df_old=pvt_old,
                df_new=pvt_new,
                title='Lab',
                path='figures/')

# include pressure
# psat = pvtc.pvt_table['p_sat']
# pvt_old.insert(0, 'p_sat', psat)
# pvt_new.insert(0, 'p_sat', psat)
#
# pvt_new.to_csv(fr'outputs/pvt_new_{type_of_data}.csv', index=False)
# pvt_old.to_csv(fr'outputs/pvt_old{type_of_data}.csv', index=False)
# pvtc.pvt_table.to_csv(fr'outputs/inputs_{type_of_data}.csv', index=False)

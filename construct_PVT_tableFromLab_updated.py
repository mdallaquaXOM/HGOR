from correlations.utils import sampling_old
from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT, excel_columns_map

# map of columns to read the Excel file
columnsNumbers, columnsNames = excel_columns_map()

# create class
pvtc = PVTCORR_HGOR(sat_pressure=None,
                    hgor=2500,
                    columns=columnsNumbers,
                    path='data',
                    files=['PVT_Data_updated'],
                    columnsNames=columnsNames,
                    categorical_columns=['well_name', 'formation', 'fluid'],
                    skiprows=0,
                    dataAugmentation=0)

# Old Properties
pvt_old = pvtc.construct_PVT_table_old()

# NEW Properties
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
    #     {'principle': 'exponential_rational_15', 'variation': 'michael'},
    'Rog':
    # {'principle': 'nasser', 'variation': 'knownPsat'},
        {'principle': 'ace', 'variation': 'ovalle'}
}
pvt_new = pvtc.construct_PVT_table_new(properties)

# Calculate metrics

# Comparing both plots.
plot_comparePVT(properties=['Rgo', 'Bo', 'visc_o', 'Rog'],
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

from correlations.utils import sampling_old
from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT

# properties and columns
knobs = dict(psat=63,
             API=62,
             gamma_s=81,  # column 81(msst?) or 100 or 99
             temperature=101,
             sep_P1=82,
             sep_T1=83,
             sep_sg1=84,
             sep_P2=85,
             sep_T2=86,
             sep_sg2=87,
             sep_P3=88,
             sep_T3=89,
             sep_sg3=90,
             sep_P4=91,
             sep_T4=92,
             sep_sg4=93,
             )

# PVTs
pvts = dict(Rgo=105,
            Rog=102,
            Bo=66,  # ???
            Bgwet=103,
            Bgdry=104,
            # Bw=,
            mu_o=72,
            # mu_g=,
            # mu_w=,
            )

extras = dict(well_name=7,
              formation=8,
              fluid=10,
              p_res=67,
              Psc=97,
              Tsc=98)

columnsNumbers = list(knobs.values()) + list(pvts.values()) + list(extras.values())
columnsNumbers = [columnsNumber - 1 for columnsNumber in columnsNumbers]

columnsNames = list(knobs.keys()) + list(pvts.keys()) + list(extras.keys())
columnsNames = [x for _, x in sorted(zip(columnsNumbers, columnsNames))]
# columnsNumbers.sort()

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
    'CGR':
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

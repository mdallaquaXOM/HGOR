from correlations.new_correlations import PVTCORR_HGOR

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat',
                             'Bo_psat', 'Bg_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data'],
                    dataAugmentation=0)

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

pvt_new = pvtc.construct_PVT_table(properties)

a = 0

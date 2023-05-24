from correlations.new_correlations import PVTCORR_HGOR

from correlations.utils import sampling, printInputValues, plot_comparePVT, concatDF
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
    #       {'principle': 'exponential_rational_15', 'variation': 'michael'},
    'muob':
    # {'principle': 'Beggs_and_Robinson', 'variation': 'original'},
        {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
    #     {'principle': 'exponential_rational_15', 'variation': 'michael'}
}

# create class
pvtc = PVTCORR_HGOR(sat_pressure=5000., Tsp=60., Psp=14.7,
                    hgor=2000.,
                    path='data',
                    files=['test_new_correlation'],
                    skiprows=[1],
                    header=0,
                    inputData=False
                    )


# Optimizer definitions
ranges = [(20, 55), (0.65, 1.2), (130, 300)]

# What was the original calculated values
pvt_old = pvtc.pvt_table
# input_original = pvtc.match_PVT_values(ranges)

# what is the new inputs values
input_star = pvtc.match_PVT_valuesHGOR(ranges,
                                       pvt_old=pvt_old,
                                       properties=properties,
                                       additional_details=True)

# Print comparison dor (API, Specific_Gravity, Temperature)
printInputValues(new=input_star)

# New Properties
pvt_new = pvtc.construct_PVT_table_new(properties, inputValues=input_star)

# Comparing both plots.
plot_comparePVT(inputs=pvtc.pvt_table[['p']],
                df_old=pvt_old,
                df_new=pvt_new,
                title='Match')

# saving the pvt table in a
pvt_new.insert(0, 'p', pvtc.pvt_table['p'])
pvt_new.to_excel(fr"outputs/macthedPVT.xlsx", index=False)


from correlations.HGOR_script import PVTCORR_HGOR

from correlations.utils import sampling_old, printInputValues, plot_comparePVT, concatDF
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
    #       {'principle': 'exponential_rational_15', 'variation': 'michael'},
    'visc_o':
    # {'principle': 'Beggs_and_Robinson', 'variation': 'original'},
        {'principle': 'Beggs_and_Robinson', 'variation': 'rs_update'},
    #     {'principle': 'exponential_rational_15', 'variation': 'michael'}
}

# create from PVT table
pvtc = PVTCORR_HGOR(sat_pressure=5000., Tsp=60., Psp=14.7,
                    hgor=2000.,
                    path='data',
                    files=['test_new_correlation'],
                    skiprows=[1],
                    header=0,
                    inputData=False
                    )
pvt_table = pvtc.pvt_table

#  Real values used to create the table
# inputs = pvt_table['p'].to_frame()
# inputs['API'] = 51.76405
# inputs['temperature'] = 261.4015
# inputs['gamma_s'] = 0.7932375
#
# pvtc_test = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
#                          hgor=2000.,
#                          data=inputs
#                          )
#
# pvt_recreated = pvtc_test.construct_PVT_table_old()

# Optimizer definitions
ranges = [(30, 55), (0.65, 1.2), (130, 300)]

# Original values calculated
input_calculated = pvtc.match_PVT_values(ranges, disp=True, maxiter=100)


# New Properties
pvt_recreated = pvtc.construct_PVT_table_new(properties, inputValues=input_calculated)


# what is the new inputs values
columnToMatch = ['Rgo', 'Bo']
input_star, _ = pvtc.match_PVT_valuesHGOR(ranges,
                                          pvt_old=pvt_table,
                                          properties=properties,
                                          additional_details=True,
                                          columnToMatch=columnToMatch,
                                          # x_start=input_calculated,
                                          disp=True)

# Print comparison dor (API, Specific_Gravity, Temperature)
printInputValues(old=input_calculated, new=input_star)

# New Properties
#todo: check this
pvt_new = pvtc.construct_PVT_table_new(properties, inputValues=input_star)

# Comparing both plots.
plot_comparePVT(inputs=pvt_table,
                df_old=pvt_recreated,
                df_new=pvt_new,
                title='Match')

# saving the pvt table in a
pvt_new.insert(0, 'psat', pvtc.pvt_table['p'])
pvt_new.to_excel(fr"outputs/macthedPVT.xlsx", index=False)

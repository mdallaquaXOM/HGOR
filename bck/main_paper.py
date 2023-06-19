import os
import matplotlib.pyplot as plt
import pandas as pd
from LGOR_script import PVTCORR_HGOR
from utils import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7, filepath=os.path.join('../..', 'Data', 'PVT_paper.xlsx'))

# # creating flag for HGOR
# plot_pairplots(pvtc.pvt_table[['p_sat', 'temperature', 'gamma_gs', 'API', 'Rs', 'HGOR']], hue='HGOR',
#                origin='paper')

# From the table
api = pvtc.pvt_table['API']
gas_gravity = pvtc.pvt_table['gamma_gs']
temperature = pvtc.pvt_table['temperature']


# Optimize Vasquez and Beggs
C_new_vasquez = pvtc.optimizeParameter(api, gas_gravity, temperature)


new_parameter = {'Vasquez_Beggs':C_new_vasquez.x}

# Calculate RS
Rs = pvtc.compute_PVT_Correlations_metrics_delete(api, gas_gravity, temperature, new_parameter)

df_Rs = pd.DataFrame(Rs)
df_Rs['HGOR'] = pvtc.pvt_table['HGOR']

plot_log_log(df_Rs, measured='Rs',
             calculated=['Vasquez_Beggs',
                         'Vasquez_Beggs_paper',
                         'Vasquez_Beggs_modified',
                         'Exponential_Rational_8',
                         'Exponential_Rational_16'],
             title='Rs (scf/stb) at saturation pressure')

plot_log_log(df_Rs, measured='Rs',
             calculated=['Vasquez_Beggs',
                         'Vasquez_Beggs_modified',
                         'Exponential_Rational_8'],
             title='Rs (scf/stb) at saturation pressure')
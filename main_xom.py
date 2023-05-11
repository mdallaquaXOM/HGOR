import os
import matplotlib.pyplot as plt
import pandas as pd
from LGOR_script import PVTCORR_HGOR
from utils import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=500, filepath=os.path.join('..', 'Data', 'PVT_Data.xlsx'))

# # creating flag for HGOR
# plot_pairplots(pvtc.pvt_table[['p_sat', 'temperature', 'gas_gravity', 'gamma_gs', 'API', 'Rs', 'HGOR']], hue='HGOR')

# From the table
api = pvtc.pvt_table['API']
gas_gravity = pvtc.pvt_table['gamma_gs']
temperature = pvtc.pvt_table['temperature']

Rs = pvtc.compute_RS_values(api, gas_gravity, temperature)

df_Rs = pd.DataFrame(Rs)
df_Rs['HGOR'] = pvtc.pvt_table['HGOR']

plot_log_log(df_Rs, measured='Rs',
             calculated=['Vasquez_Beggs',
                         'Vasquez_Beggs_modified',
                         'Exponential_Rational_8',
                         'Exponential_Rational_16'],
             title='Rs (scf/stb) at saturation pressure')

a = 0

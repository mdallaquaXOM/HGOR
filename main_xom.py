from LGOR_script import PVTCORR_HGOR
from utils import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat', 'Bo_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper'])


# Perform EDA
# EDA_seaborn(pvtc.pvt_table)
# EDA_plotly(pvtc.pvt_table)


# Optimize Vasquez and Beggs
C_new_vasquez = pvtc.optimizeParameter(metric_func='LSE')

new_parameter = {'Vasquez_Beggs': C_new_vasquez.x}

# Calculate RS
Rs = pvtc.compute_RS_values(new_parameter)

df_Rs = pd.DataFrame(Rs)
df_Rs['HGOR'] = pvtc.pvt_table['HGOR']

plot_log_log(df_Rs, measured='Rs',
             calculated=['VB_original',
                         'VB_modified',
                         'VB_paper',
                         'Exp_Rational_8',
                         'Exp_Rational_16_paper',
                         'Exp_Rational_16_new'],
             title='Rs (scf/stb) at saturation pressure')

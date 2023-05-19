import pickle
from correlations.new_correlations import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter
from correlations.definitions import *

pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                    hgor=2000,
                    columns=['temperature', 'API', 'gamma_s', 'Rs', 'p_sat', 'Bo_psat', 'Bg_psat', 'visc_o_psat'],
                    path='data',
                    files=['PVT_Data', 'PVT_paper']
                    )

source_plots = 'PVT_Data'

# Calculate RS
pvt_df, pvt_metrics = pvtc.compute_PVT_values(source=source_plots)

# plots
for index, row in pvt_df.iterrows():

    df = row.to_frame().T
    df = df.explode(['actual', 'calculated']).reset_index(drop=True)
    df_hgor = pvtc.pvt_table[pvtc.pvt_table['source'] == source_plots]
    df['HGOR'] = df_hgor['HGOR']

    plot_log_log(df, measured='actual',
                 calculated=['calculated'],
                 metrics_df=pvt_metrics.loc[index].to_frame().T,
                 title=index, property=index)

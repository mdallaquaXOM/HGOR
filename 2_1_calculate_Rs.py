from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import *

source_curve = 'PVT_Data_updated'
source_curve = 'PVT_Data'

# create class
if source_curve == 'PVT_Data_updated':
    # map of columns to read the Excel file
    columnsNumbers, columnsNames = excel_columns_map()

    pvtc = PVTCORR_HGOR(sat_pressure=None,
                        hgor=2500,
                        columns=columnsNumbers,
                        path='data',
                        files=[source_curve],
                        columnsNames=columnsNames,
                        categorical_columns=['well_name', 'formation', 'fluid'],
                        skiprows=0,
                        dataAugmentation=0)
else:
    pvtc = PVTCORR_HGOR(sat_pressure=None, Tsp=60, Psp=14.7,
                        hgor=2500,
                        columns=['temperature', 'API', 'gamma_s', 'Rgo', 'psat', 'Bo', 'visc_o'],
                        path='data',
                        files=['PVT_Data', 'PVT_paper'],
                        dataAugmentation=0)

# properties to plot
properties = {'Rgo': [
    {'principle': 'vasquez_beggs', 'variation': 'original'},
    {'principle': 'ace', 'variation': 'mine'},
    {'principle': 'datadriven', 'variation': 'ann'},
    {'principle': 'datadriven', 'variation': 'randomforest'},
    # {'principle': 'vasquez_beggs', 'variation': 'optimized'},
    # {'principle': 'vasquez_beggs', 'variation': 'meija'},
    {'principle': 'exponential_rational_8', 'variation': 'blasingame'},
    {'principle': 'exponential_rational_8', 'variation': 'optimized'},
    # {'principle': 'exponential_rational_16', 'variation': 'blasingame'},
    # {'principle': 'exponential_rational_16', 'variation': 'michael'},
],
}

# Calculate RS
# pvt_prop, pvt_metrics = pvtc.compute_PVT_Correlations_metrics_delete(properties,
#                                                                      new_parameters=new_parameters,
#                                                                      source=source_curve,
#                                                                         )
pvt_prop = pvtc.compute_PVT_Correlations(properties,
                                         source=source_curve)

# calculate metrics
pvt_metrics, pvt_prop = calculateMetrics(pvtc.pvt_table, pvt_prop, source=source_curve)

# plot
colums2plot = pvt_prop['Rgo'].drop(['measured', 'HGOR'], axis=1).columns.values

plot_properties(pvt_prop['Rgo'], measured='measured',
                calculated=colums2plot,
                metrics_df=pvt_metrics['Rgo'],
                title='Rs (scf/stb) at saturation pressure')

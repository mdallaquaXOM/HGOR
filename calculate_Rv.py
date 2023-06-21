from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import *

# from page 108

# range of variables
p_bounds = [4000., 7500]

bounds = {'temperature': [75, 225],
          'gamma_s': [0.65, 0.85],
          'API': [40, 60],
          }

# New Correlations
properties = {
    'Rog': [
        {'principle': 'nasser', 'variation': 'VO_unknownPsat'},
        {'principle': 'nasser', 'variation': 'VO_knownPsat'},
        {'principle': 'nasser', 'variation': 'GC_unknownPsat'},
        {'principle': 'nasser', 'variation': 'GC_knownPsat'},
        {'principle': 'ace', 'variation': 'ovalle'},
    ],
}

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

# NEW Properties
pvt_prop = pvtc.compute_PVT_Correlations(properties)

# calculate metrics
pvt_metrics, pvt_prop = calculateMetrics(pvtc.pvt_table, pvt_prop)

# plot
colums2plot = pvt_prop['Rog'].drop(['measured', 'HGOR'], axis=1).columns.values

plot_properties(pvt_prop['Rog'], measured='measured',
                calculated=colums2plot,
                metrics_df=pvt_metrics['Rog'],
                property='Rv',
                title='Rv (stb/MMscf) at saturation pressure')

# input_table = pvt_prop.pvt_table

# for property_, correlations in pvt_prop.items():
#     plot_synthetic_data(correlations, input_table, name=property_)
#
#     plot_synthetic_data(correlations, input_table, name=property_, hueplot='sample')

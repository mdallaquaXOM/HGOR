from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import EDA_seaborn, EDA_plotly, plot_comparePVT, excel_columns_map

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
#
EDA_seaborn(pvtc.pvt_table, hue='HGOR')
EDA_seaborn(pvtc.pvt_table, hue='fluid')

a = 0
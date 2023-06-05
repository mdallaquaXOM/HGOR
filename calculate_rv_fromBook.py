# Replicate Example 4 from El-Banbi, Ahmed, Ahmed Alzahabi, and Ahmed El-Maraghi. 2018. PVT Property Correlations -
# Selection and Estimation. Saint Louis: Elsevier

from correlations.HGOR_script import PVTCORR_HGOR
from correlations.utils import *
from correlations.optimizer import optimizeParameter
from correlations.definitions import *

# from page 108
Rvi = 230  # STB/MMscf
gamma_sp_1 = 0.786
gamma_sp_2 = 0.788
Tr = 323


# New Correlations  - select JUST ONE !!!!!!
properties = {
    # 'Rs':
    #     {'principle': 'nasser', 'variation': 'GC'},
    # {'principle': 'nasser', 'variation': 'VO'},
    # 'Bo':
    #     {'principle': 'nasser', 'variation': 'GC'},
    # {'principle': 'nasser', 'variation': 'VO'},
    'Rv':
        {'principle': 'nasser', 'variation': 'GC'},
    # {'principle': 'nasser', 'variation': 'VO'},
    # 'Bg':
    #     {'principle': 'nasser', 'variation': 'GC'},
    # {'principle': 'nasser', 'variation': 'VO'},
}

pvtc = PVTCORR_HGOR(sat_pressure=None,
                    Tsp1=110, Psp1=800,
                    Tsp2=75, Psp2=60,
                    hgor=2000,
                    path='data',
                    files=['PVT_Book'],
                    correct_specific_gravity=False)

pvt_new = pvtc.construct_PVT_table_new(properties)

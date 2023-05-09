import os
import matplotlib.pyplot as plt
import pandas as pd
from LGOR_script import PVTCORR_HGOR


if __name__ == '__main__':
    sat_pressure = 5123.2

    pvtc = PVTCORR_HGOR(sat_pressure=sat_pressure, Tsp=60, Psp=500, filepath=os.path.join('..', 'Data', 'PVT Data.xlsx'))

    #          API     gas_gravity  temperature
    ranges = [(50, 55), (0.85, 1), (130, 230)]
    api, gas_gravity, temperature = pvtc.match_PVT_values(ranges, additional_details=True)
    dct = pvtc.compute_PVT_values(api, gas_gravity, temperature)
    # dct = pvtc.compute_PVT_values(48.35,   0.817, 136.31)
    # dct = pvtc.compute_PVT_values(65.01, 1.0648, 237.28)
    df = pd.DataFrame(dct)
    df[df['pressure'] <= sat_pressure].to_csv('plu_PVT_500.csv')
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_Rgo', 'Calculated_Rgo'], ylim=[0, 6000])
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_Bo', 'Calculated_Bo'], ylim=[0.5, 3.5])
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_Bg', 'Calculated_Bg'], ylim=[0, 0.05])
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_Bw', 'Calculated_Bw'], ylim=[0.5, 1.5])
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_vo', 'Calculated_vo'], ylim=[0, 5])
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_vg', 'Calculated_vg'], ylim=[0, 0.1])
    df[df['pressure'] <= sat_pressure].plot(x='pressure', y=['Actual_vw', 'Calculated_vw'], ylim=[0, 1])
plt.show()
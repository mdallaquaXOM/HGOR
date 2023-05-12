from .utils import metrics
import numpy as np
import os
import pandas as pd
from .LGOR_script import PVTCORR


class PVTCORR_HGOR(PVTCORR):
    def __init__(self, path, files, columns=None, hgor=2000, **kwargs):

        super().__init__(**kwargs)

        pvt_tables = []
        for file in files:
            filepath = os.path.join(path, file + '.xlsx')

            pvt_table_i = pd.read_excel(filepath, header=1, usecols=columns)
            pvt_table_i['source'] = file
            pvt_tables.append(pvt_table_i)

        pvt_table = pd.concat(pvt_tables)

        # Calculate corrected gas gravity (?)
        if 'gamma_c' not in pvt_table.columns:
            api = pvt_table['API']
            gas_gravity = pvt_table['gamma_s']

            gamma_gs = self._computeGasGravityAtSeparatorConditions(gas_gravity, api)
            pvt_table['gamma_c'] = gamma_gs

        # Assign flag for HGOR
        pvt_table['HGOR'] = False
        pvt_table.loc[pvt_table['Rs'] > hgor, 'HGOR'] = True

        self.pvt_table = pvt_table

    def _computeSolutionGasOilRatio(self, api, temperature,
                                    pressure, gas_gravity, method=None,
                                    coefficients=None):
        # standard values
        if method is None:
            method = {'principle': 'vasquez_beggs', 'variation': 'original'}
        if coefficients is None:
            coefficients = [0.0178, 1.187, 23.931, 47.]

        if method['principle'] == "vasquez_beggs":
            if method['variation'] == "original":
                conditions = [api <= 30, (api > 30)]

                C1_choices = [0.0362, 0.0178]
                C2_choices = [1.0937, 1.187]
                C3_choices = [25.724, 23.931]

            elif method['variation'] == 'extra':

                conditions = [api <= 30, (api > 30) & (api < coefficients[3]), api >= coefficients[3]]

                C1_choices = [0.0362, 0.0178, coefficients[0]]
                C2_choices = [1.0937, 1.187, coefficients[1]]
                C3_choices = [25.724, 23.931, coefficients[2]]

            else:
                raise ValueError(f"Unknown {method['variation']} for calculating Rs ")

            C1 = np.select(conditions, C1_choices, default=np.nan)
            C2 = np.select(conditions, C2_choices, default=np.nan)
            C3 = np.select(conditions, C3_choices, default=np.nan)

            a = pressure ** C2
            b = np.exp((C3 * api) / (temperature + 459.67))

            Rs = C1 * a * gas_gravity * b

        elif method['principle'] == "vasquez_beggs_paper":
            C1 = np.where(api <= 30, 1.091e+5, 3.405e+6)
            C2 = np.where(api <= 30, 2.3913, 2.7754)
            C3 = np.where(api <= 30, 1.e-6, 1.e-6)

            a = pressure ** C3
            b = - (C2 * api) / (temperature + 459.67)

            Rs = (gas_gravity * a * (10 ** b)) / C1

        elif method['principle'] == "exponential_rational_8":
            C = [9.021, -0.119, 2.221, -.531, .144, -1.842e-2, 12.802, 8.309]

            a = C[0] + C[1] * np.log(temperature)
            b = C[2] + C[3] * np.log(api)
            d = C[6] + C[7] * np.log(gas_gravity)

            K = (a / np.log(pressure) - 1.) * (b * d) ** (-1)

            ln_Rs = (K - C[4]) / C[5]

            Rs = np.exp(ln_Rs)

        elif method['principle'] == "exponential_rational_16":
            if method['variation'] == 'new_paper':
                # New paper
                C = [
                    7.258546e-1, -4.562008e-2,
                    3.198814e00, -3.994698e-1,
                    -1.483415e-1, 3.550853e-1,
                    2.914460e00, 4.402225e-1,
                    -1.791551e-1, 6.955443e-1,
                    -8.172007e-1, 4.229810e-1,
                    -5.612631e-1, 4.735904e-02,
                    4.746990e-02, -2.515009e-01
                ]
            else:  ## Blasingane paper
                C = [0.858, -7.881e-2,
                     3.198, -.457,
                     .146, .322,
                     3.172, 1.015,
                     -.34, .54,
                     -.665, .458,
                     -.545, 3.343e-2,
                     .454, -.281
                     ]

            ln_t = np.log(temperature)
            ln_api = np.log(api)
            ln_gas_gravity = np.log(gas_gravity)

            a = C[0] + C[1] * ln_t
            e = C[8] + C[9] * ln_t

            b = C[2] + C[3] * ln_api
            f = C[10] + C[11] * ln_api

            c = C[6] + C[7] * ln_gas_gravity
            g = C[14] + C[15] * ln_gas_gravity

            A = a * b * c
            B = e * f * g

            K = np.log(pressure) * B / A

            ln_Rs = (K * C[12] - C[4]) / (C[5] - K * C[13])

            Rs = np.exp(ln_Rs)

        else:
            raise ValueError(f'Unknown method ({method}) for calculating Rs ')

        return Rs

    def compute_RS_values(self, new_parameter=None, source=None):

        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask]

        # recover inputs
        api = df['API']
        gas_gravity = df['gamma_c']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])

        rs = np.array(df['Rs'])

        # New correlations
        rs_vb_orig = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                      method={'principle': 'vasquez_beggs', 'variation': 'original'})

        rs_vb_mod = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                     method={'principle': 'vasquez_beggs', 'variation': 'extra'},
                                                     coefficients=new_parameter['Vasquez_Beggs'])

        rs_vb_paper = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                       method={'principle': 'vasquez_beggs_paper'})

        rs_exp_rat_8 = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                        method={'principle': 'exponential_rational_8'})

        rs_exp_rat_16_paper = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                               method={'principle': 'exponential_rational_16',
                                                                       'variation': 'blasingame'})

        rs_exp_rat_16_new = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                             method={'principle': 'exponential_rational_16',
                                                                     'variation': 'new_paper'})

        comparison_dict = {'pressure': p_sat,
                           'temperature': temperature,
                           'gas_gravity': gas_gravity,
                           'api': api,
                           'Rs': rs,
                           'VB_original': rs_vb_orig,
                           'VB_modified': rs_vb_mod,
                           'VB_paper': rs_vb_paper,
                           'Exp_Rational_8': rs_exp_rat_8,
                           'Exp_Rational_16_paper': rs_exp_rat_16_paper,
                           'Exp_Rational_16_new': rs_exp_rat_16_new}

        # Old correlations
        # comparison_dict = {'Vasquez_Beggs': []}
        # for p_i, api_i, gas_gravity_i, temperature_i in zip(p_sat, api, gas_gravity, temperature):
        #     self.sat_pressure = p_i
        #     rs_vb = super()._computeSolutionGasOilRatio(api_i, temperature_i, p_i, gas_gravity_i)
        #     comparison_dict['Vasquez_Beggs'].append(rs_vb)

        metrics_ = {'VB_original': metrics(rs, rs_vb_orig),
                    'VB_modified': metrics(rs, rs_vb_mod),
                    'VB_paper': metrics(rs, rs_vb_paper),
                    'Exp_Rational_8': metrics(rs, rs_exp_rat_8),
                    'Exp_Rational_16_paper': metrics(rs, rs_exp_rat_16_paper),
                    'Exp_Rational_16_new': metrics(rs, rs_exp_rat_16_new)}

        # treat outputs
        comparison_df = pd.DataFrame.from_dict(comparison_dict)
        comparison_df['HGOR'] = df['HGOR']

        metrics_df = pd.DataFrame.from_dict(metrics_, orient='index')
        metrics_df = metrics_df.round(2)

        return comparison_df, metrics_df



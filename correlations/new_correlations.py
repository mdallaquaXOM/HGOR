from .utils import metrics
import numpy as np
import os
import pandas as pd
from .LGOR_script import PVTCORR
from .definitions import *
import pickle


class PVTCORR_HGOR(PVTCORR):
    def __init__(self, path, files, columns=None, hgor=2000, dataAugmentation=None, **kwargs):

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

        if dataAugmentation is not None:
            hgor = pvt_table[pvt_table['HGOR'] == True]
            for i in range(dataAugmentation):
                pvt_table = pd.concat([pvt_table, hgor])

        self.pvt_table = pvt_table

    def _computeSolutionGasOilRatio(self, api, temperature,
                                    pressure, gas_gravity, method=None,
                                    parameters=None):
        # standard values
        if method is None:
            method = {'principle': 'vasquez_beggs', 'variation': 'original'}

        principle = method['principle'].lower()
        variation = method['variation'].lower()

        if principle == "vasquez_beggs":
            if variation == "original":
                conditions = [api <= 30, (api > 30)]

                C1_choices = C_vasquez_beggs['C1']
                C2_choices = C_vasquez_beggs['C2']
                C3_choices = C_vasquez_beggs['C3']

            elif variation == 'optimized':

                conditions = [api <= 30,
                              (api > 30) & (api < parameters[3]),
                              api >= parameters[3]]

                C1_choices = [*C_vasquez_beggs['C1'], parameters[0]]
                C2_choices = [*C_vasquez_beggs['C2'], parameters[1]]
                C3_choices = [*C_vasquez_beggs['C3'], parameters[2]]

            elif variation == 'meija':
                C1 = np.where(api <= 30, *C_vasquez_beggs_meija_martinez['C1'], )
                C2 = np.where(api <= 30, *C_vasquez_beggs_meija_martinez['C2'], )
                C3 = np.where(api <= 30, *C_vasquez_beggs_meija_martinez['C3'], )

                a = pressure ** C3
                b = - (C2 * api) / (temperature + 459.67)

                Rs = (gas_gravity * a * (10 ** b)) / C1

                return Rs

            else:
                msg = f"Unknown variation:{variation} for method:{principle} calculating Rs"
                raise ValueError(msg)

            C1 = np.select(conditions, C1_choices, default=np.nan)
            C2 = np.select(conditions, C2_choices, default=np.nan)
            C3 = np.select(conditions, C3_choices, default=np.nan)

            a = pressure ** C2
            b = np.exp((C3 * api) / (temperature + 459.67))

            Rs = C1 * a * gas_gravity * b

        elif principle == "exponential_rational_8":
            if variation == 'optimized':
                C = parameters
            elif variation == 'blasingame':
                C = C_exp_rat_8_blasingame
            elif variation == 'meija':
                C = C_exp_rat_16_meija_martinez
            else:  # blasingame paper
                raise Exception(f'No variation: {variation} for method: {principle}')

            a = C[0] + C[1] * np.log(temperature)
            b = C[2] + C[3] * np.log(api)
            d = C[6] + C[7] * np.log(gas_gravity)

            K = (a / np.log(pressure) - 1.) * ((b * d) ** (-1))

            ln_Rs = (K - C[4]) / C[5]

            Rs = np.exp(ln_Rs)

        elif principle == "exponential_rational_16":
            if variation == 'michael':
                # New paper
                C = C_exp_rat_16_michael
            elif variation == 'blasingame':  # Blasingame paper
                C = C_exp_rat_16_blasingame
            else:  # method['variation'] == 'optimized':
                C = parameters

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

            ln_pb = np.log(pressure)

            ln_Rs = (A * C[4] - ln_pb * (1 + B * C[12])) / (-A * C[5] + B * C[13] * ln_pb)

            Rs = np.exp(ln_Rs)

        else:
            raise ValueError(f'Unknown method ({method}) for calculating Rs ')

        return Rs

    def compute_RS_values(self, correlations=None, new_parameter=None, source=None):

        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask].reset_index(drop=True)

        if correlations is None:
            correlations = [{'principle': 'vasquez_beggs', 'variation': 'original'}]

        # recover inputs
        api = df['API']
        gas_gravity = df['gamma_c']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])
        rs = np.array(df['Rs'])

        # New correlations
        rs_values = {'method': ['measured'], 'values': [rs]}
        metrics_dic = {'method': [], 'values': []}
        for correlation in correlations:
            rs_temp = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                       method=correlation,
                                                       parameters=new_parameter[correlation['principle']])

            method = correlation['principle'] + '_' + correlation['variation']
            rs_values['method'].append(method)
            rs_values['values'].append(rs_temp)

            metrics_dic['method'].append(method)
            metrics_dic['values'].append(metrics(rs, rs_temp))

        # convert things to dataframe
        comparison_df = pd.DataFrame(rs_values['values'], index=rs_values['method']).T
        metrics_df = pd.DataFrame(metrics_dic['values'], index=metrics_dic['method'])
        metrics_df = metrics_df.round(2)

        # add HGOR flag
        comparison_df['HGOR'] = df['HGOR']

        return comparison_df, metrics_df

    def _compute_bubblePressure(self, api, temperature, rs, gas_gravity
                                , method=None, parameters=None):
        if method == 'vasquez_beggs':
            conditions = [api <= 30, (api > 30)]

            C1_choices = C_vasquez_beggs['C1']
            C2_choices = C_vasquez_beggs['C2']
            C3_choices = C_vasquez_beggs['C3']

            C1 = np.select(conditions, C1_choices, default=np.nan)
            C2 = np.select(conditions, C2_choices, default=np.nan)
            C3 = np.select(conditions, C3_choices, default=np.nan)

            a = np.exp(C3 * api / (temperature + 459.67))
            b = C1 * gas_gravity * a

            pb = (rs / b) ** (1 / C2)

        elif method == 'vasquez_beggs_blasingame':
            conditions = [api <= 30, (api > 30)]

            C1_choices = C_vasquez_beggs_meija_martinez['C1']
            C2_choices = C_vasquez_beggs_meija_martinez['C2']
            C3_choices = C_vasquez_beggs_meija_martinez['C3']

            C1 = np.select(conditions, C1_choices, default=np.nan)
            C2 = np.select(conditions, C2_choices, default=np.nan)
            C3 = np.select(conditions, C3_choices, default=np.nan)

            a = 10. ** (C2 * api / (temperature + 459.67))
            b = C1 * (rs / gas_gravity) * a
            pb = b ** (C3)

        elif method == 'exponential_rational_8':
            C = C_exp_rat_8_blasingame

            ln_t = np.log(temperature)
            ln_api = np.log(api)
            ln_gas_gravity = np.log(gas_gravity)
            ln_gas_rs = np.log(rs)

            a = C[0] + C[1] * ln_t
            b = C[2] + C[3] * ln_api
            c = C[4] + C[5] * ln_gas_rs
            d = C[6] + C[7] * ln_gas_gravity

            ln_pb = a / (1 + b * c * d)

            pb = np.exp(ln_pb)

        elif method[:23] == 'exponential_rational_16':
            method_ = method.split('_')[-1]
            if method_ == 'blasingame':
                C = C_exp_rat_16_blasingame
            else:
                C = parameters

            ln_t = np.log(temperature)
            ln_api = np.log(api)
            ln_gas_gravity = np.log(gas_gravity)
            ln_gas_rs = np.log(rs)

            a = C[0] + C[1] * ln_t
            e = C[8] + C[9] * ln_t

            b = C[2] + C[3] * ln_api
            f = C[10] + C[11] * ln_api

            c = C[4] + C[5] * ln_gas_rs
            g = C[12] + C[13] * ln_gas_rs

            d = C[6] + C[7] * ln_gas_gravity
            h = C[14] + C[15] * ln_gas_gravity

            ln_pb = (a * b * c * d) / (1. + e * f * g * h)

            pb = np.exp(ln_pb)
        else:
            raise Exception(f'wrong method: {method}')

        return pb

    def compute_pb_values(self, source=None):
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
        pb_vb_orig = self._compute_bubblePressure(api, temperature, rs, gas_gravity,
                                                  method='vasquez_beggs')

        pb_vb_paper = self._compute_bubblePressure(api, temperature, rs, gas_gravity,
                                                   method='vasquez_beggs_blasingame')

        pb_exp_rat_8_paper = self._compute_bubblePressure(api, temperature, rs, gas_gravity,
                                                          method='exponential_rational_8')

        pb_exp_rat_16_paper = self._compute_bubblePressure(api, temperature, rs, gas_gravity,
                                                           method='exponential_rational_16_blasingame')

        pb_exp_rat_16_ed = self._compute_bubblePressure(api, temperature, rs, gas_gravity,
                                                        method='exponential_rational_16_opt',
                                                        parameters=C_exp_rat_16_michael)

        comparison_dict = {'pressure': p_sat,
                           'temperature': temperature,
                           'gas_gravity': gas_gravity,
                           'api': api,
                           'Rs': rs,
                           'VB_original': pb_vb_orig,
                           'VB_paper': pb_vb_paper,
                           'Exp_Rational_8_paper': pb_exp_rat_8_paper,
                           'Exp_Rational_16_paper': pb_exp_rat_16_paper,
                           # 'Exp_Rational_16_opt': pb_exp_rat_16_opt,
                           'Exp_Rational_16_Ed': pb_exp_rat_16_ed
                           }

        metrics_ = {'VB_original': metrics(p_sat, pb_vb_orig),
                    'VB_paper': metrics(p_sat, pb_vb_paper),
                    'Exp_Rational_8_paper': metrics(p_sat, pb_exp_rat_8_paper),
                    'Exp_Rational_16_paper': metrics(p_sat, pb_exp_rat_16_paper),
                    'Exp_Rational_16_Ed': metrics(p_sat, pb_exp_rat_16_ed),
                    # 'Exp_Rational_16_optimized': metrics(rs, rs_exp_rat_16_opt)
                    }

        # treat outputs
        comparison_df = pd.DataFrame.from_dict(comparison_dict)
        comparison_df['HGOR'] = df['HGOR']

        metrics_df = pd.DataFrame.from_dict(metrics_)
        metrics_df = metrics_df.round(2)

        return comparison_df, metrics_df

    def compute_PVT_values(self, source=None, correlations=None):
        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask]

        # recover inputs
        api = df['API']
        gas_gravity_s = df['gamma_s']
        gas_gravity = df['gamma_c']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])

        # recover outputs
        rgo = np.array(df['Rs'])
        bo = np.array(df['Bo_psat'])
        bg = np.array(df['Bg_psat'])
        Visc_o = np.array(df['visc_o_psat'])

        # optimized parameter
        new_parameter = pickle.load(open(r"optimizedParam/opt_results.pickle", "rb"))

        # Gor at saturation pressure - Rs
        rgo_c = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                 method={'principle': 'exponential_rational_8',
                                                         'variation': 'optimize'},
                                                 parameters=new_parameter['exponential_rational_8'])

        # old correlations
        bo_c = self._computeLiveOilFVF(api, temperature, p_sat, gas_gravity, rgo_c)

        bo_c_old = []
        bg_c = []
        bw_c = []
        Visc_o_c = []
        Visc_g_c = []
        Visc_w_c = []

        # Old Correlations Below
        for p_sat_i, temperature_i, gas_gravity_i, api_i \
                in zip(p_sat, temperature, gas_gravity_s, api):
            self.sat_pressure = p_sat_i

            bo_c_old.append(super()._computeLiveOilFVF(api_i, temperature_i, p_sat_i, gas_gravity_i))
            bg_c.append(self.computeDryGasFVF(p_sat_i, temperature_i, gas_gravity_i))
            bw_c.append(self.computeWaterFVF(temperature_i, p_sat_i))
            Visc_o_c.append(self.computeLiveOilViscosity(api_i, temperature_i, p_sat_i, gas_gravity_i))
            Visc_g_c.append(self.computeDryGasViscosity(temperature_i, p_sat_i, gas_gravity_i))
            Visc_w_c.append(self.computerWaterViscosity(p_sat_i, temperature_i))

        bo_c_old = np.asarray(bo_c_old, dtype=np.float64)
        bg_c = np.asarray(bg_c, dtype=np.float64)
        bw_c = np.asarray(bw_c, dtype=np.float64)
        Visc_o_c = np.asarray(Visc_o_c, dtype=np.float64)
        Visc_g_c = np.asarray(Visc_g_c, dtype=np.float64)
        Visc_w_c = np.asarray(Visc_w_c, dtype=np.float64)

        comparison_dict = {'Rgo': {'actual': rgo, 'calculated': rgo_c},
                           'Bo': {'actual': bo, 'calculated': bo_c},
                           'Bg': {'actual': bg, 'calculated': bg_c},
                           'Visc_o': {'actual': Visc_o, 'calculated': Visc_o_c},
                           'bw': {'actual': None, 'calculated': bw_c},
                           'Visc_g': {'actual': None, 'calculated': Visc_g_c},
                           'Visc_w': {'actual': None, 'calculated': Visc_w_c},
                           }

        metrics_ = {'Rgo': metrics(rgo, rgo_c),
                    'Bo': metrics(bo, bo_c),
                    'Bg': metrics(bg, bg_c),
                    'Visc_o': metrics(Visc_o, Visc_o_c),
                    }

        # treat outputs
        comparison_df = pd.DataFrame.from_dict(comparison_dict, orient='index')
        # comparison_df['HGOR'] = df['HGOR']

        metrics_df = pd.DataFrame.from_dict(metrics_, orient='index')
        metrics_df = metrics_df.round(2)

        return comparison_df, metrics_df

    def _computeLiveOilFVF(self, api, temperature, pressure, gas_gravity,
                           rs, method=None, parameters=None):
        # Vazquez and Beggs, 1980 (Default in EMPower)

        # standard values
        if method is None:
            method = {'principle': 'vasquez_beggs', 'variation': 'original'}

        if method['principle'] == "vasquez_beggs":
            if method['variation'] == "blasingame":
                C4_choices = C_vasquez_beggs_meija_martinez['C4']
                C5_choices = C_vasquez_beggs_meija_martinez['C5']
                C6_choices = C_vasquez_beggs_meija_martinez['C6']
            else:
                C4_choices = C_vasquez_beggs['C4']
                C5_choices = C_vasquez_beggs['C5']
                C6_choices = C_vasquez_beggs['C6']

            conditions = [api <= 30, (api > 30)]

            C4 = np.select(conditions, C4_choices, default=np.nan)
            C5 = np.select(conditions, C5_choices, default=np.nan)
            C6 = np.select(conditions, C6_choices, default=np.nan)

            Tr = temperature - 60.
            a = Tr * (api / gas_gravity)

            Bo = 1. + C4 * rs + (C5 + C6 * rs) * a

        else:
            raise Exception(f'Method {method} not implemented')

        return Bo

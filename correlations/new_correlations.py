from .utils import metrics
import numpy as np
import os
import pandas as pd
from .LGOR_script import PVTCORR
from .definitions import *
import pickle


class PVTCORR_HGOR(PVTCORR):
    def __init__(self, path='', files='', data=None, columns=None, hgor=2000, dataAugmentation=None, **kwargs):

        super().__init__(**kwargs)

        if data is None:
            pvt_tables = []
            for file in files:
                filepath = os.path.join(path, file + '.xlsx')

                pvt_table_i = pd.read_excel(filepath, header=1, usecols=columns)
                pvt_table_i['source'] = file
                pvt_tables.append(pvt_table_i)

            pvt_table = pd.concat(pvt_tables)
        else:
            pvt_table = data

        # Calculate corrected gas gravity (?)
        if 'gamma_c' not in pvt_table.columns:
            api = pvt_table['API']
            gas_gravity = pvt_table['gamma_s']

            gamma_gs = self._computeGasGravityAtSeparatorConditions(gas_gravity, api)
            pvt_table['gamma_c'] = gamma_gs

        # Assign flag for HGOR
        if 'Rs' in pvt_table:
            pvt_table['HGOR'] = False
            pvt_table.loc[pvt_table['Rs'] > hgor, 'HGOR'] = True

            if dataAugmentation is not None:
                hgor = pvt_table[pvt_table['HGOR'] == True]
                for i in range(dataAugmentation):
                    pvt_table = pd.concat([pvt_table, hgor])

        self.pvt_table = pvt_table

    @staticmethod
    def _computeRsAtSatPress(api, temperature,
                             pressure, gas_gravity, method=None,
                             parameters=None):
        # standard values
        if method is None:
            method = {'principle': 'exponential_rational_8', 'variation': 'optimized'}

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
                if parameters is None:
                    new_parameter = pickle.load(open(r"optimizedParam/opt_results.pickle", "rb"))
                    C = new_parameter['Rs']['exponential_rational_8']
                else:
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

    @staticmethod
    def _compute_bubblePressure(api, temperature, rs, gas_gravity, method=None, parameters=None):
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

    def _computeBoAtSatPres(self, api, temperature, psat, gas_gravity, method=None, rs=None, rs_best_correlation=None):

        # standard values
        if method is None:
            method = {'principle': 'vasquez_beggs', 'variation': 'original'}

        principle = method['principle'].lower()
        variation = method['variation'].lower()

        if rs is None:
            if principle == "exponential_rational_15":
                rs = self._computeRsAtSatPress(api, temperature, psat, gas_gravity, method=rs_best_correlation)
            elif principle == "vasquez_beggs":
                if variation == "rs_update":
                    rs = self._computeRsAtSatPress(api, temperature, psat, gas_gravity, method=rs_best_correlation)
                else:
                    rs = self._computeRsAtSatPress(api, temperature, psat, gas_gravity,
                                                   method={'principle': 'vasquez_beggs', 'variation': 'original'})

        if principle == "vasquez_beggs":
            if variation == "meija":
                C4_choices = C_vasquez_beggs_meija_martinez['C4']
                C5_choices = C_vasquez_beggs_meija_martinez['C5']
                C6_choices = C_vasquez_beggs_meija_martinez['C6']
            elif variation == "original" or variation == "rs_update":
                C4_choices = C_vasquez_beggs['C4']
                C5_choices = C_vasquez_beggs['C5']
                C6_choices = C_vasquez_beggs['C6']
            else:
                raise Exception(f'Bob method {principle} with variation ({variation}) not coded')

            conditions = [api <= 30, (api > 30)]

            C4 = np.select(conditions, C4_choices, default=np.nan)
            C5 = np.select(conditions, C5_choices, default=np.nan)
            C6 = np.select(conditions, C6_choices, default=np.nan)

            Tr = temperature - 60.
            a = Tr * (api / gas_gravity)

            Bob = 1. + C4 * rs + (C5 + C6 * rs) * a

        elif principle == "exponential_rational_15":
            C = C_exp_rat_15_Bob

            ln_t = np.log(temperature)
            ln_api = np.log(api)
            ln_rs = np.log(rs)
            ln_gas_gravity = np.log(gas_gravity)
            ln_psat = np.log(psat)

            ln_t_2 = ln_t ** 2
            ln_api_2 = ln_api ** 2
            ln_rs_2 = ln_rs ** 2
            ln_gas_gravity_2 = ln_gas_gravity ** 2
            ln_psat_2 = ln_psat ** 2

            ln_Bob = (C[0] + C[1] * ln_t + C[2] * ln_t_2) * \
                     (C[3] + C[4] * ln_api + C[5] * ln_api_2) * \
                     (C[6] + C[7] * ln_rs + C[8] * ln_rs_2) * \
                     (C[9] + C[10] * ln_gas_gravity + C[11] * ln_gas_gravity_2) * \
                     (C[12] + C[13] * ln_psat + C[14] * ln_psat_2)

            Bob = np.exp(ln_Bob)

        else:
            raise Exception(f'Method {method} not implemented')

        return Bob

    def _computeMuobAtSatPres(self, api, temperature, psat, gas_gravity, method=None, Rso=None,
                              rs_best_correlation=None):
        # standard values
        if method is None:
            method = {'principle': 'Beggs_and_Robinson'}

        principle = method['principle'].lower()
        variation = method['variation'].lower()

        if Rso is None:
            if principle == "exponential_rational_15":
                Rso = self._computeRsAtSatPress(api, temperature, psat, gas_gravity, method=rs_best_correlation)
            elif principle == "beggs_and_robinson":
                if variation == "rs_update":
                    Rso = self._computeRsAtSatPress(api, temperature, psat, gas_gravity, method=rs_best_correlation)
                else:
                    Rso = self._computeRsAtSatPress(api, temperature, psat, gas_gravity,
                                                    method={'principle': 'vasquez_beggs', 'variation': 'original'})

        if principle == "beggs_and_robinson":
            # Beggs and Robinson, 1975 Defailt EMPower
            Visc_oil = self.computeDeadOilViscosity(api, temperature)
            a = 10.715 * ((Rso + 100.0) ** (-0.515))
            b = 5.44 * ((Rso + 150.0) ** (-0.338))
            Visc_oil = a * (Visc_oil ** b)

        elif principle == "exponential_rational_15":
            C = C_exp_rat_15_muob

            ln_t = np.log(temperature)
            ln_api = np.log(api)
            ln_rs = np.log(Rso)
            ln_gas_gravity = np.log(gas_gravity)
            ln_psat = np.log(psat)

            ln_t_2 = ln_t ** 2
            ln_api_2 = ln_api ** 2
            ln_rs_2 = ln_rs ** 2
            ln_gas_gravity_2 = ln_gas_gravity ** 2
            ln_psat_2 = ln_psat ** 2

            ln_muob = (C[0] + C[1] * ln_t + C[2] * ln_t_2) * \
                      (C[3] + C[4] * ln_api + C[5] * ln_api_2) * \
                      (C[6] + C[7] * ln_rs + C[8] * ln_rs_2) * \
                      (C[9] + C[10] * ln_gas_gravity + C[11] * ln_gas_gravity_2) * \
                      (C[12] + C[13] * ln_psat + C[14] * ln_psat_2)

            Visc_oil = np.exp(ln_muob)

        else:
            raise Exception(f'Method {method} not implemented')

        return Visc_oil

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

        metrics_dict = {'VB_original': metrics(p_sat, pb_vb_orig),
                        'VB_paper': metrics(p_sat, pb_vb_paper),
                        'Exp_Rational_8_paper': metrics(p_sat, pb_exp_rat_8_paper),
                        'Exp_Rational_16_paper': metrics(p_sat, pb_exp_rat_16_paper),
                        'Exp_Rational_16_Ed': metrics(p_sat, pb_exp_rat_16_ed),
                        # 'Exp_Rational_16_optimized': metrics(rs, rs_exp_rat_16_opt)
                        }

        # treat outputs
        comparison_df = pd.DataFrame.from_dict(comparison_dict)
        comparison_df['HGOR'] = df['HGOR']

        metrics_df = pd.DataFrame.from_dict(metrics_dict).T
        metrics_df = metrics_df.round(2)

        return comparison_df, metrics_df

    def compute_PVT_Correlations(self, properties, rs_best_correlation=None,
                                 new_parameters=None, source=None):

        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask].reset_index(drop=True)

        # recover inputs
        api = df['API']
        gas_gravity = df['gamma_c']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])
        rs = np.array(df['Rs'])
        Bo_psat = np.array(df['Bo_psat'])
        visc_o_psat = np.array(df['visc_o_psat'])

        # New correlations
        comparison_star = {}
        metrics_star = {}

        new_parameter = None
        newparameter_prop_correl = None

        for property_, correlations in properties.items():
            # get the new parameters for the property in question
            if new_parameters is not None:
                if property_ in new_parameters:
                    new_parameter = new_parameters[property_]

            # get measured value
            if property_ == 'Rs':
                value_measured = rs
            elif property_ == 'Bo':
                value_measured = Bo_psat
            elif property_ == 'muob':
                value_measured = visc_o_psat
            else:
                raise Exception(f'Not able to get the measured values')

            prop_values = {'method': ['measured'], 'values': [value_measured]}
            metrics_dic = {'method': [], 'values': []}

            for correlation in correlations:
                if new_parameter is not None:
                    newparameter_prop_correl = new_parameter[correlation['principle']]

                # function to call will depend on the property_
                if property_ == 'Rs':
                    temp = self._computeRsAtSatPress(api, temperature, p_sat, gas_gravity,
                                                     method=correlation,
                                                     parameters=newparameter_prop_correl)
                elif property_ == 'Bo':
                    temp = self._computeBoAtSatPres(api, temperature, p_sat, gas_gravity,
                                                    method=correlation,
                                                    rs_best_correlation=rs_best_correlation)

                elif property_ == 'muob':
                    temp = self._computeMuobAtSatPres(api, temperature, p_sat, gas_gravity,
                                                      method=correlation,
                                                      rs_best_correlation=rs_best_correlation)
                else:
                    raise Exception(f'Property unknown {property_}')

                method = correlation['principle'] + '_' + correlation['variation']
                prop_values['method'].append(method)
                prop_values['values'].append(temp)

                metrics_dic['method'].append(method)
                metrics_dic['values'].append(metrics(value_measured, temp))

            # convert things to dataframe
            comparison_df = pd.DataFrame(prop_values['values'], index=prop_values['method']).T
            metrics_df = pd.DataFrame(metrics_dic['values'], index=metrics_dic['method'])
            metrics_df = metrics_df.round(2)

            # add HGOR flag
            comparison_df['HGOR'] = df['HGOR']

            comparison_star[property_] = comparison_df
            metrics_star[property_] = metrics_df

        return comparison_star, metrics_star

    def compute_PVT_Correlations_v2(self, properties, rs_best_correlation=None,
                                    new_parameters=None, source=None):

        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask].reset_index(drop=True)

        # recover inputs
        api = df['API']
        gas_gravity = df['gamma_c']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])

        # New correlations
        comparison_star = {}

        new_parameter = None
        newparameter_prop_correl = None

        for property_, correlations in properties.items():
            # get the new parameters for the property in question
            if new_parameters is not None:
                if property_ in new_parameters:
                    new_parameter = new_parameters[property_]

            prop_values = {'method': [], 'values': []}

            for correlation in correlations:
                if new_parameter is not None:
                    newparameter_prop_correl = new_parameter[correlation['principle']]

                # function to call will depend on the property_
                if property_ == 'Rs':
                    temp = self._computeRsAtSatPress(api, temperature, p_sat, gas_gravity,
                                                     method=correlation,
                                                     parameters=newparameter_prop_correl)
                elif property_ == 'Bo':
                    temp = self._computeBoAtSatPres(api, temperature, p_sat, gas_gravity,
                                                    method=correlation,
                                                    rs_best_correlation=rs_best_correlation)

                elif property_ == 'muob':
                    temp = self._computeMuobAtSatPres(api, temperature, p_sat, gas_gravity,
                                                      method=correlation,
                                                      rs_best_correlation=rs_best_correlation)
                else:
                    raise Exception(f'Property unknown {property_}')

                method = correlation['principle'] + '_' + correlation['variation']
                prop_values['method'].append(method)
                prop_values['values'].append(temp)

            # convert things to dataframe
            comparison_df = pd.DataFrame(prop_values['values'], index=prop_values['method']).T

            # # add HGOR flag
            # comparison_df['HGOR'] = df['HGOR']

            comparison_star[property_] = comparison_df

        return comparison_star

    def construct_PVT_table(self, properties, rs_best_correlation=None,
                            new_parameters=None, source=None):

        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask].reset_index(drop=True)

        # recover inputs
        api = df['API']
        gas_gravity_s = df['gamma_s']
        gas_gravity = df['gamma_c']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])

        # Pvt properties
        rgo_c = self._computeRsAtSatPress(api, temperature, p_sat, gas_gravity,
                                          method=properties['Rs'])

        bo_c = self._computeBoAtSatPres(api, temperature, p_sat, gas_gravity,
                                        method=properties['Bo'],
                                        rs=rgo_c)

        Visc_o_c = self._computeMuobAtSatPres(api, temperature, p_sat, gas_gravity,
                                              method=properties['muob'],
                                              Rso=rgo_c)

        # Old Correlations Below
        # rgo_c = []
        # bo_c = []
        bg_c = []
        bw_c = []
        # Visc_o_c = []
        Visc_g_c = []
        Visc_w_c = []

        for p_sat_i, temperature_i, gas_gravity_i, api_i \
                in zip(p_sat, temperature, gas_gravity_s, api):
            self.sat_pressure = p_sat_i

            # rgo_c.append(super()._computeSolutionGasOilRatio(api_i, temperature_i, p_sat_i, gas_gravity_i))
            # bo_c.append(super()._computeLiveOilFVF(api_i, temperature_i, p_sat_i, gas_gravity_i))
            bg_c.append(super().computeDryGasFVF(p_sat_i, temperature_i, gas_gravity_i))
            bw_c.append(super().computeWaterFVF(temperature_i, p_sat_i))
            # Visc_o_c.append(super().computeLiveOilViscosity(api_i, temperature_i, p_sat_i, gas_gravity_i))
            Visc_g_c.append(super().computeDryGasViscosity(temperature_i, p_sat_i, gas_gravity_i))
            Visc_w_c.append(super().computerWaterViscosity(p_sat_i, temperature_i))

        pvt_dic = {
            'Rgo': rgo_c,
            'Bo': bo_c,
            'Bg': bg_c,
            'Bw': bw_c,
            'Visc_o': Visc_o_c,
            'Visc_g': Visc_g_c,
            'Visc_w': Visc_w_c,
        }

        # treat outputs
        pvt_df = pd.DataFrame.from_dict(pvt_dic)

        return pvt_df

    def construct_PVT_table_old(self, source=None):
        df = self.pvt_table

        # filter by source
        if source is not None:
            mask = df['source'] == source
            df = df[mask]

        # recover inputs
        api = df['API']
        gas_gravity_s = df['gamma_s']
        temperature = df['temperature']
        p_sat = np.array(df['p_sat'])

        # recover outputs
        # rgo = np.array(df['Rs'])
        # bo = np.array(df['Bo_psat'])
        # bg = np.array(df['Bg_psat'])
        # Visc_o = np.array(df['visc_o_psat'])

        rgo_c = []
        bo_c = []
        bg_c = []
        bw_c = []
        Visc_o_c = []
        Visc_g_c = []
        Visc_w_c = []

        # Old Correlations Below
        for p_sat_i, temperature_i, gas_gravity_i, api_i \
                in zip(p_sat, temperature, gas_gravity_s, api):
            self.sat_pressure = p_sat_i

            rgo_c.append(super()._computeSolutionGasOilRatio(api_i, temperature_i, p_sat_i, gas_gravity_i))
            bo_c.append(super()._computeLiveOilFVF(api_i, temperature_i, p_sat_i, gas_gravity_i))
            bg_c.append(super().computeDryGasFVF(p_sat_i, temperature_i, gas_gravity_i))
            bw_c.append(super().computeWaterFVF(temperature_i, p_sat_i))
            Visc_o_c.append(super().computeLiveOilViscosity(api_i, temperature_i, p_sat_i, gas_gravity_i))
            Visc_g_c.append(super().computeDryGasViscosity(temperature_i, p_sat_i, gas_gravity_i))
            Visc_w_c.append(super().computerWaterViscosity(p_sat_i, temperature_i))

        bo_c = np.asarray(bo_c, dtype=np.float64)
        bg_c = np.asarray(bg_c, dtype=np.float64)
        bw_c = np.asarray(bw_c, dtype=np.float64)
        Visc_o_c = np.asarray(Visc_o_c, dtype=np.float64)
        Visc_g_c = np.asarray(Visc_g_c, dtype=np.float64)
        Visc_w_c = np.asarray(Visc_w_c, dtype=np.float64)

        pvt_dic = {
            'Rgo': rgo_c,
            'Bo': bo_c,
            'Bg': bg_c,
            'Bw': bw_c,
            'Visc_o': Visc_o_c,
            'Visc_g': Visc_g_c,
            'Visc_w': Visc_w_c,
        }

        # treat outputs
        pvt_df = pd.DataFrame.from_dict(pvt_dic)
        # comparison_df['HGOR'] = df['HGOR']

        return pvt_df

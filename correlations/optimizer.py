import numpy as np
import scipy as sp
from .utils import metrics


def optimizeParameter(self, source=None, metric_func='LSE'):
    df = self.pvt_table

    # filter by soruce
    if source is not None:
        mask = df['source'] == source
        df = df[mask]

    # recover inputs
    api = df['API']
    gas_gravity = df['gamma_c']
    temperature = df['temperature']
    p_sat = np.array(df['p_sat'])

    # For metric evaluation
    rs_measured = np.array(df['Rs'])

    # objective function
    def obj_function(coefficients):
        rs_calc = self._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                   method={'principle': 'vasquez_beggs', 'variation': 'extra'},
                                                   coefficients=coefficients)
        metrics_ = metrics(rs_measured, rs_calc)
        obj = metrics_[metric_func]

        return obj

    #                       C1         C2        C3         C5
    range_of_values = [(1e-2, 8e-2), (1., 2.), (15., 30), (40., 50)]
    x0 = [0.0178, 1.187, 23.931, 47.]
    C_new = sp.optimize.differential_evolution(
        obj_function, range_of_values,
        seed=100,
        x0=x0,
        strategy='best2exp')

    print(C_new)

    return C_new

import numpy as np
import scipy as sp
from .utils import metrics
from scipy.optimize import minimize, Bounds, differential_evolution


def optimizeParameter(pvt_class, source=None, metric_func='LSE', type=3):
    # type
    # 1 - Trust-Region Constrained Algorithm (method='trust-constr')
    # 2 - Sequential Least SQuares Programming (SLSQP) Algorithm (method='SLSQP')
    # 3 - Global optimization

    df = pvt_class.pvt_table

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
        rs_calc = pvt_class._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                        method={'principle': 'vasquez_beggs', 'variation': 'extra'},
                                                        coefficients=coefficients)
        metrics_ = metrics(rs_measured, rs_calc)
        obj = metrics_[metric_func]

        return obj

    #                    lower                 upper
    #                 C1   C2   C3   C4     C1   C2  C3  C4
    bounds = Bounds([1e-2, 1., 15., 40.], [8e-2, 2., 30, 50])

    x_start = np.array([0.0178, 1.187, 23.931, 47.])

    # Optimization
    if type == 1 | type == 2:
        if type == 1:
            print('Optimizer: Trust-Region Constrained Algorithm')
            method = 'trust-constr'
        else:
            print('Optimizer: Sequential Least SQuares Programming (SLSQP) Algorithm')
            method = 'SLSQP'

        x_new = minimize(obj_function, x_start, method=method, bounds=bounds, tol=1e-8,
                         options={'disp': True})
    else:
        print('Optimizer: Global optimization')
        x_new = differential_evolution(obj_function, bounds=bounds, x0=x_start, tol=1e-8)

    print(f'Best set of parameters: {x_new.x}')

    return x_new

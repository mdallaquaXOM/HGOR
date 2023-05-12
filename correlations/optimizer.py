import numpy as np
import scipy as sp
from .utils import metrics
from scipy.optimize import minimize, Bounds, differential_evolution


def optimizeParameter(pvt_class,
                      source=None,
                      metric_func='LSE',
                      algorithm=3,
                      correlation_method=None,
                      bounds=None,
                      x_start=None):
    # type
    # 1 - Global optimization
    # 2 - Trust-Region Constrained Algorithm (method='trust-constr')
    # 3 - Sequential Least SQuares Programming (SLSQP) Algorithm (method='SLSQP')

    df = pvt_class.pvt_table

    # filter by source
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
    def obj_function(parameters):
        rs_calc = pvt_class._computeSolutionGasOilRatio(api, temperature, p_sat, gas_gravity,
                                                        method=correlation_method,
                                                        parameters=parameters)
        metrics_ = metrics(rs_measured, rs_calc)
        obj = metrics_[metric_func]

        return obj

    # Optimization
    if algorithm == 1:
        print('\tOptimizer: Global optimization')
        x_new = differential_evolution(obj_function, bounds=Bounds(*bounds), x0=x_start, tol=1e-8)
    else:
        if algorithm == 2:
            print('\tOptimizer: Trust-Region Constrained Algorithm')
            method = 'trust-constr'
        elif algorithm == 3:
            print('\tOptimizer: Sequential Least SQuares Programming (SLSQP) Algorithm')
            method = 'SLSQP'
        elif algorithm == 4:
            print("\tUnconstrained minimization: Broyden-Fletcher-Goldfarb-Shanno algorithm (method='BFGS')")
            method = 'BFGS'
        else:
            method = None

        if bounds is not None:
            bounds = Bounds(*bounds)

        x_new = minimize(obj_function, x_start, method=method, tol=1e-8, bounds=bounds,
                         options={'disp': True})

    print(f'\tBest set of parameters: {x_new.x}')

    return x_new

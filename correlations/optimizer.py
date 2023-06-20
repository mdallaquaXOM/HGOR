import numpy as np
import scipy as sp
from .utils import metrics
from scipy.optimize import minimize, Bounds, differential_evolution


def optimizeParameter(pvt_class,
                      opt_equation='Rs',
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

    if bounds is not None:
        bounds = Bounds(*bounds)

    # recover inputs
    api = df['API']
    gas_gravity = df['gamma_c']
    temperature = df['temperature']
    p_sat = np.array(df['psat'])

    # For metric evaluation
    rs_measured = np.array(df['Rgo'])

    if opt_equation == 'Rs':
        # objective function
        def obj_function(parameters):
            rs_calc = pvt_class._computeRsAtSatPress(api, temperature, p_sat, gas_gravity,
                                                     method=correlation_method,
                                                     parameters=parameters)
            metrics_ = metrics(rs_measured, rs_calc)
            obj = metrics_[metric_func]

            return obj
    elif opt_equation == 'pb':
        def obj_function(parameters):
            pb_calc = pvt_class._compute_bubblePressure(api, temperature, rs_measured, gas_gravity,
                                                        method=correlation_method,
                                                        parameters=parameters)
            metrics_ = metrics(p_sat, pb_calc)
            obj = metrics_[metric_func]
            # if obj == np.inf:
            #     obj = pvt_class.LARGE

            return obj

    # Optimization
    if algorithm == 1:
        print('\tOptimizer: Global optimization')
        x_new = differential_evolution(obj_function, bounds=bounds, x0=x_start, tol=1e-8, strategy='best2exp')
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

        x_new = minimize(obj_function, x_start, method=method, tol=1e-8, bounds=bounds,
                         options={'disp': True})

    print(f'\tBest set of parameters: {x_new.x}')

    return x_new


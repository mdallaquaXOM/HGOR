import numpy as np
import scipy as sp
from .utils import metrics
from scipy.optimize import minimize, Bounds, differential_evolution


def optimizeParameter(pvt_class,
                      opt_equation='Rgo',
                      source=None,
                      metric_func='LSE',
                      algorithm=3,
                      correlation_method=None,
                      bounds=None,
                      x_start=None,
                      correct_gamma=True,
                      verbose=False,
                      callback=None):
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
    temperature = df['temperature']
    p_sat = np.array(df['psat'])

    if correct_gamma:
        gas_gravity = df['gamma_c']
    else:
        gas_gravity = df['gamma_s']

    # For metric evaluation
    measured = np.array(df[opt_equation])

    # objective function
    def obj_function(parameters, info=None):
        if opt_equation == 'Rgo':
            calc = pvt_class._computeRsAtSatPress(api, temperature, p_sat, gas_gravity,
                                                  method=correlation_method,
                                                  parameters=parameters)
        elif opt_equation == 'pb':
            rs_measured = np.array(df['Rgo'])
            calc = pvt_class._compute_bubblePressure(api, temperature, rs_measured, gas_gravity,
                                                     method=correlation_method,
                                                     parameters=parameters)
        elif opt_equation == 'Rog':
            separators = {'specific_gravity': [df['sep_sg1'], df['sep_sg2'], df['sep_sg3'], df['sep_sg4']],
                          'pressure': [df['sep_P1'], df['sep_P2'], df['sep_P3'], df['sep_P4']],
                          'temperature': [df['sep_T1'], df['sep_T2'], df['sep_T3'], df['sep_T4']],
                          }
            calc = pvt_class._computeVaporizedOilGas(API=api, temperature=temperature, pressure=p_sat,
                                                     gas_gravity=gas_gravity, method=correlation_method,
                                                     separators=separators,
                                                     parameters=parameters
                                                     )

        else:
            raise Exception(f'opt_equation not coder {opt_equation}')

        metrics_ = metrics(measured, calc)
        obj = metrics_[metric_func]

        if info is not None:
            if info['Nfeval'] % 1000 == 0:
                strFormat = '{0:4d}   {1: 3.4e}   {2: 3.4e}   {3: 3.4e}   {4: 3.4e}   {5: 3.4e}   {6: 3.4e}   {7: ' \
                            '.4e}   ' \
                            '{8: 3.4e}   {9: 3.4e}'
                print(strFormat.format(info['Nfeval'], *parameters, obj))
            info['Nfeval'] += 1

        return obj

    # Optimization
    if algorithm == 1:

        if callback:
            print('\tOptimizer: Global optimization')
            print()
            strFormat = '{0:4s}   {1:11s}   {2:11s}   {3:11s}   {4:11s}   {5:11s}   {6:11s}   {7:11s}   {8:11s}   {' \
                        '9:11s}'
            print(strFormat.format('Iter', ' C0', ' C1', ' C2', ' C3', ' C4', ' C5', ' C6', ' C7', ' f(X)'))

            args = ({'Nfeval': 0},)
        else:
            args = None

        x_new = differential_evolution(obj_function, bounds, x0=x_start, tol=1e-8, strategy='best2exp',
                                       disp=verbose, args=args)
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
                         options={'disp': verbose})

    print(f'\tBest set of parameters: {x_new.x}')

    return x_new

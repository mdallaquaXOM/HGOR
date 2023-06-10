from ace import ace
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re


def read_column_data_from_xlxs(fname, columns=None, independent_var=None, log=None, header=1, skiprows=None):
    pvt_table = pd.read_excel(fname, header=header, usecols=columns, skiprows=skiprows)

    # convert to log if asked
    if log is not None:
        for var_log in log:
            pvt_table[var_log] = np.log(pvt_table[var_log])
            new_column_name = f'{var_log}_ln'
            pvt_table = pvt_table.rename(columns={var_log: new_column_name})

            if independent_var == var_log:
                independent_var = new_column_name

    y_values = pvt_table[[independent_var]]
    x_values = pvt_table.drop([independent_var], axis=1)

    var_names = {'independent': x_values.columns.to_list(),
                 'dependent': y_values.columns.to_list()}

    y_values = y_values.values.flatten()
    x_values = x_values.values.T.tolist()

    return x_values, y_values, var_names


def read_column_data_from_dat(fname, independent_var=None, log=None):
    """
    Read data from a simple text file.

    Format should be just numbers.
    First column is the dependent variable. others are independent.
    Whitespace delimited.

    Returns
    -------
    x_values : list
        List of x columns
    y_values : list
        list of y values

    """
    datafile = open(fname)
    datarows = []

    for iline, line in enumerate(datafile):
        if iline == 0:
            headers = line.split()
        else:
            datarows.append([float(li) for li in line.split()])
    datacols = list(zip(*datarows))

    # convert to log if asked
    if log is not None:
        for var_log in log:
            log_var_idx = headers.index(var_log)
            log_var = np.log(datacols[log_var_idx])
            datacols[log_var_idx] = log_var
            # datacols[log_var_idx] = tuple(map(tuple, log_var))

            headers[log_var_idx] = f'ln_{headers[log_var_idx]}'

            if independent_var == var_log:
                independent_var = headers[log_var_idx]

    y_var = headers.index(independent_var)

    y_values = datacols.pop(y_var)
    x_values = datacols

    y_name = headers.pop(y_var)

    var_names = {'independent': headers,
                 'dependent': y_name}

    return x_values, y_values, var_names


def postProcessing(ace_model, order_indep=2, order_dep=2, var_names=None, fname=None):
    print()
    print('Coefficients for regression:')

    # Regression for independent variables
    str_names = []
    beta_independet = []

    nVariables = len(ace_model.x)
    if not isinstance(order_indep, list):
        order_indep = [order_indep] * nVariables

    for i in range(nVariables):
        x = np.asarray(ace_model.x[i]).reshape(-1, 1)
        y = ace_model.x_transforms[i].reshape(-1, 1)

        poly_reg = PolynomialFeatures(degree=order_indep[i], include_bias=True)
        X_poly = poly_reg.fit_transform(X=x)
        # pol_reg = LinearRegression()
        # pol_reg.fit(X_poly, y)
        # a = pol_reg.intercept_
        # beta = pol_reg.coef_

        beta = np.linalg.solve(X_poly.T @ X_poly, X_poly.T @ y)
        # print(f"\t {var_names['independent'][i]}_Tr = {['' for i, beta_i in enumerate(beta) ]}")
        # print(f"\t {var_names['independent'][i]}_Tr = {beta.reshape(1, -1)}")
        var_name = var_names['independent'][i]
        str_names.append(f"\t {var_name}_Tr = "
                         f""
                         f""
                         f"{' + '.join(map(str, [f'{beta_i[0]:.4e}*{var_name}**{j}' for j, beta_i in enumerate(beta)]))}")
        print(str_names[i])
        beta_independet.append(beta)

    # Regression for dependent variables
    x = np.asarray(ace_model.y_transform).reshape(-1, 1)
    y = ace_model.y.reshape(-1, 1)

    poly_reg = PolynomialFeatures(degree=order_dep, include_bias=True)
    X_poly = poly_reg.fit_transform(X=x)

    beta_dep = np.linalg.solve(X_poly.T @ X_poly, X_poly.T @ y)
    var_name = var_names['dependent']
    str_names.append('')
    str_names.append(f"\t {var_name[0]} = "f""
                     f"{' + '.join(map(str, [f'{beta_i[0]:.4e}*Sum_Tr**{j}' for j, beta_i in enumerate(beta_dep)]))}")
    print(str_names[-1])

    coefficients = {'independent': beta_independet,
                    'dependent': beta_dep}

    # str_names = [re.sub('\t ', '', str_name) for str_name in str_names]
    full_var = []
    for list_var in var_names.values():
        for var_name in list_var:
            full_var.append(var_name)

    order = order_indep.copy()
    order.append(order_dep)

    # save coefficient to txt
    betas = beta_independet
    betas.append(beta_dep)
    with open(fname, 'w') as f:
        for i_beta, beta in enumerate(betas):
            line = [f'{betai[0]:.4e}' for betai in beta]
            line.append('\n')
            f.write(f'{full_var[i_beta]}, {order[i_beta]}: ' + ', '.join(line))
            # f.write(', '.join(beta))
    # ace_model._write_columns(fname, beta_independet, [beta_dep])

    # np.savetxt(fname, betas, fmt='%.4e', delimiter=', ')

    return coefficients


def plot_input(ace_model, fname='ace_input.png', var_names=None):
    """Plot the transforms."""
    if not plt:
        raise ImportError('Cannot plot without the matplotlib package')
    plt.rcParams.update({'font.size': 8})
    plt.figure()
    num_cols = len(ace_model.x) // 2 + 1
    for i in range(len(ace_model.x)):
        plt.subplot(num_cols, 2, i + 1)
        plt.plot(ace_model.x[i], ace_model.y, '.')
        if var_names is not None:
            plt.xlabel(f"{var_names['independent'][i]}")
            plt.ylabel(f"{var_names['dependent']}")
        else:
            plt.xlabel('x{0}'.format(i))
            plt.ylabel('y')

    plt.tight_layout()

    if fname:
        plt.savefig(fname)
    else:
        plt.show()


def plot_optimalRegression(ace_model, fname='ace_transforms.png', var_names=None):
    """Plot the transforms."""
    if not plt:
        raise ImportError('Cannot plot without the matplotlib package')
    plt.rcParams.update({'font.size': 8})
    plt.figure()

    sumTR = 0
    for i in range(len(ace_model.x)):
        sumTR += ace_model.x_transforms[i]

    plt.plot(sumTR, ace_model.y_transform, '.')
    if var_names is not None:
        plt.xlabel('Sum_Tr_Independent')
        plt.ylabel(f"{var_names['dependent']}_Tr")
    else:
        plt.xlabel('Sum_Phi')
        plt.ylabel('theta')

    # Regression for transformed variables
    x = sumTR.reshape(-1, 1)
    y = ace_model.y_transform

    poly_reg = PolynomialFeatures(degree=1, include_bias=True)
    X_poly = poly_reg.fit_transform(X=x)
    beta_transformed = np.linalg.solve(X_poly.T @ X_poly, X_poly.T @ y)
    print()
    print(f'Transformed correlation: Theta = {beta_transformed[1]:.4e} Phi + {beta_transformed[0]:.4e}')

    # Inserting 45 degree
    x_min = np.amin(sumTR)
    x_max = np.amax(sumTR)

    x_45 = np.linspace(x_min, x_max)
    plt.plot(x_45, x_45, '--', label='45 degrees')

    plt.tight_layout()

    if fname:
        plt.savefig(fname)
        return None
    return plt


def plot_transforms(ace_model, fname='ace_transforms.png', var_names=None, regression_coeff=None):
    """Plot the transforms."""
    if not plt:
        raise ImportError('Cannot plot without the matplotlib package')
    plt.rcParams.update({'font.size': 8})
    plt.figure()
    num_cols = len(ace_model.x) // 2 + 1
    for i in range(len(ace_model.x)):
        plt.subplot(num_cols, 2, i + 1)
        plt.plot(ace_model.x[i], ace_model.x_transforms[i], '.', label='Phi {0}'.format(i))

        if regression_coeff is not None:
            coefficient = regression_coeff['independent'][i]
            minx = np.amin(ace_model.x[i])
            maxx = np.amax(ace_model.x[i])

            x_var = np.linspace(start=minx, stop=maxx).reshape(-1, 1)

            order = len(coefficient) - 1
            poly_reg = PolynomialFeatures(degree=order, include_bias=True)
            X_poly = poly_reg.fit_transform(X=x_var)

            y_pred = X_poly @ coefficient

            plt.plot(x_var, y_pred, label='regression')

        if var_names is not None:
            plt.xlabel(f"{var_names['independent'][i]}")
            plt.ylabel(f"{var_names['independent'][i]}_Tr")
        else:
            plt.xlabel('x{0}'.format(i))
            plt.ylabel('phi{0}'.format(i))
        plt.legend()

    # plotting output
    plt.subplot(num_cols, 2, i + 2)  # pylint: disable=undefined-loop-variable
    plt.plot(ace_model.y_transform, ace_model.y, '.', label='Theta')

    if regression_coeff is not None:
        coefficient = regression_coeff['dependent']
        minx = np.amin(ace_model.y_transform)
        maxx = np.amax(ace_model.y_transform)

        x_var = np.linspace(start=minx, stop=maxx).reshape(-1, 1)

        order = len(coefficient) - 1
        poly_reg = PolynomialFeatures(degree=order, include_bias=True)
        X_poly = poly_reg.fit_transform(X=x_var)

        y_pred = X_poly @ coefficient

        plt.plot(x_var, y_pred, label='regression')

    if var_names is not None:
        plt.ylabel(f"{var_names['dependent']}")
        plt.xlabel(f"{var_names['dependent']}_Tr")
    else:
        plt.ylabel('y')
        plt.xlabel('theta')
    plt.tight_layout()
    plt.legend()

    if fname:
        plt.savefig(fname)
        return None
    return plt


def plot_optimalInvTranform(ace_model, fname='ace_transforms.png', var_names=None):
    """Plot the transforms."""
    if not plt:
        raise ImportError('Cannot plot without the matplotlib package')
    plt.rcParams.update({'font.size': 8})
    plt.figure()

    plt.plot(ace_model.y_transform, ace_model.y, '.')

    if var_names is not None:
        plt.xlabel(f"{var_names['dependent']}_Tr")
        plt.ylabel(f"{var_names['dependent']}")
    else:
        plt.xlabel('theta')
        plt.ylabel('y')
    plt.tight_layout()

    if fname:
        plt.savefig(fname)
        return None
    return plt

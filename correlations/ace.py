from ace import ace
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import numpy as np
import matplotlib.pyplot as plt


def postProcessing(ace_model, order_indep=2, order_dep=2, var_names=None, fname=None):
    print('Coefficients for regression:')
    # Regression for independent variables
    str_names = []
    coeff_ind = []
    for i in range(len(ace_model.x)):
        x = np.asarray(ace_model.x[i]).reshape(-1, 1)
        y = ace_model.x_transforms[i].reshape(-1, 1)

        poly_reg = PolynomialFeatures(degree=order_indep, include_bias=True)
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
                         f"{' + '.join(map(str, [f'{beta_i[0]:.4e}*{var_name}^{j}' for j, beta_i in enumerate(beta)]))}")
        print(str_names[i])
        coeff_ind.append(beta)

    # Regression for dependent variables
    x = np.asarray(ace_model.y_transform).reshape(-1, 1)
    y = ace_model.y.reshape(-1, 1)

    poly_reg = PolynomialFeatures(degree=order_dep, include_bias=True)
    X_poly = poly_reg.fit_transform(X=x)
    # pol_reg = LinearRegression()
    # pol_reg.fit(X_poly, y)
    # a = pol_reg.intercept_
    # beta = pol_reg.coef_

    beta = np.linalg.solve(X_poly.T @ X_poly, X_poly.T @ y)
    var_name = var_names['dependent']
    str_names.append('')
    str_names.append(f"\t {var_name} = "f""
                     f"{' + '.join(map(str, [f'{beta_i[0]:.4e}*Sum_Tr^{j}' for j, beta_i in enumerate(beta)]))}")
    print(str_names[-1])

    # print(f"\t {var_names['dependent']} = {beta.reshape(1, -1)}")

    coeff_dep = beta

    coefficients = {'independent': coeff_ind,
                    'dependent': coeff_dep}

    # save coefficient to txt
    if fname is not None:
        with open(fname, 'w') as output_file:
            for line in str_names:
                output_file.write(f"{line}\n")

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

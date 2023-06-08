import matplotlib.pyplot as plt
import numpy as np


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




from ace import model
from correlations.ace import postProcessing, plot_transforms, plot_optimalRegression, read_column_data_from_xlxs

# get x and y


myace = model.Model()

x, y, var_names = read_column_data_from_xlxs(fname=r'data/PVT_Data.xlsx',
                                             columns=['temperature', 'API', 'gamma_s', 'Rgo', 'p'],
                                             independent_var='Rgo',
                                             log=['API', 'gamma_s', 'Rgo', 'p'])

myace.build_model_from_xy(x, y)

myace.ace.write_transforms_to_file(fname=r'outputs/ace/dg.txt')

# post-processing
coefficients = postProcessing(myace.ace, order_indep=[2, 4, 2, 2], order_dep=2, var_names=var_names,
                              fname=r'outputs/ace/coefficients.txt')

# plotting
# ace.plot_input(myace.ace, fname=r'outputs/dg_inputs.pdf', var_names=var_names)
plot_transforms(myace.ace, fname=r'outputs/ace/rs_transforms.png', var_names=var_names, regression_coeff=coefficients)
plot_optimalRegression(myace.ace, fname=r'outputs/ace/rs_optimalTrans.png', var_names=var_names)

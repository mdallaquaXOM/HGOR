from ace import model
from correlations.ace import postProcessing, plot_transforms, plot_optimalRegression
from correlations.utils import read_column_data_from_dat

# get x and y


myace = model.Model()

x, y, var_names = read_column_data_from_dat(fname=r'Ace/data/pvt.dat', independent_var='pb', log=['pb'])
myace.build_model_from_xy(x, y)

myace.ace.write_transforms_to_file(fname=r'outputs/ace/dg.txt')

# post-processing
coefficients = postProcessing(myace.ace, order_indep=4, order_dep=2, var_names=var_names,
                              fname=r'outputs/ace/coefficients.txt')

# plotting
# ace.plot_input(myace.ace, fname=r'outputs/dg_inputs.pdf', var_names=var_names)
plot_transforms(myace.ace, fname=r'outputs/ace/dg_transforms.pdf', var_names=var_names, regression_coeff=coefficients)
plot_optimalRegression(myace.ace, fname=r'outputs/ace/dg_transformations.pdf', var_names=var_names)

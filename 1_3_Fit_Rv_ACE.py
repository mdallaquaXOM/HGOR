from ace import model
from correlations.ace import postProcessing, plot_transforms, plot_optimalRegression, read_column_data_from_xlxs
from correlations.utils import excel_columns_map

# get data
columnsNumbers, columnsNames = excel_columns_map()
x, y, vars_ = read_column_data_from_xlxs(fname=r'data/PVT_Data_updated.xlsx',
                                         columns=columnsNumbers,
                                         columnsNames=columnsNames,
                                         independent_var='Rog',
                                         dependent_vars=['API', 'gamma_s', 'temperature', 'psat'],
                                         # log=['API', 'gamma_s', 'temperature',  'psat','Rog'],
                                         log=['psat', 'Rog'],
                                         )

# create instance
myace = model.Model()
myace.build_model_from_xy(x, y)

myace.ace.write_transforms_to_file(fname=r'outputs/ace/rv_transformed_variables.txt')
plot_transforms(myace.ace, fname=r'outputs/ace/rv_transforms.png', variables=vars_)
plot_optimalRegression(myace.ace, fname=r'outputs/ace/rv_optimalTrans.png', var_names=vars_['var_names'])

# post-processing
coefficients = postProcessing(myace.ace, order_indep=[2, 2, 5, 2], order_dep=1, variables=vars_,
                              fname=r'outputs/ace/rv_coefficients.txt')

# plotting
# ace.plot_input(myace.ace, fname=r'outputs/dg_inputs.pdf', var_names=var_names)
plot_transforms(myace.ace, fname=r'outputs/ace/rv_transforms_fitted.png', variables=vars_,
                regression_coeff=coefficients)

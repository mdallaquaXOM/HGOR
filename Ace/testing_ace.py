from ace.samples import wang04
from ace import model, ace

# To use, get some sample data:
x, y = wang04.build_sample_ace_problem_wang04(N=200)




# and run:
myace = model.Model()
myace.build_model_from_xy(x, y)
# myace.build_model_from_txt(fname = 'myinput.txt')
myace.eval([0.1, 0.2, 0.5, 0.3, 0.5])

# plotting
ace.plot_transforms(myace.ace, fname=r'outputs/tutorial.pdf')
myace.ace.write_transforms_to_file(fname=r'outputs/tutorial.txt')

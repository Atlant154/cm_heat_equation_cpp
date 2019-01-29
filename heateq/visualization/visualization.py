from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import json
import os
import fnmatch

viz_files = []

for filename in os.listdir("."):
    if fnmatch.fnmatch(filename, 'visualization_*.txt'):
        viz_files.append(filename)

if not viz_files:
    print("No file to visualize!")
    exit(1)

for file in viz_files:

    print(file)

    with open(file) as handle:
        json_data = json.loads(handle.read())

    x_left_bound  = np.float64(json_data['x_left_bound'])
    x_right_bound = np.float64(json_data['x_right_bound'])
    h_num         = np.float64(json_data['h_number'])
    x_lin         = np.linspace(x_left_bound, x_right_bound, h_num)

    t_left_bound  = np.float64(json_data['t_left_bound'])
    t_right_bound = np.float64(json_data['t_right_bound'])
    tau_num       = np.float64(json_data['tau_number'])
    y_lin         = np.float64(np.linspace(t_left_bound, t_right_bound, tau_num))

    X, Y = np.meshgrid(x_lin, y_lin)
    Z = np.array(json_data['result'])

    print(X.shape)
    print(Y.shape)
    print(Z.shape)

    figure = plt.figure()
    axis = figure.gca(projection='3d')
    surface = axis.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=1, antialiased=False)

    axis.zaxis.set_major_locator(LinearLocator(10))
    axis.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    figure.colorbar(surface, shrink=0.5, aspect=5)

    axis.set_xlabel('X')
    axis.set_ylabel('Time')
    axis.set_zlabel(json_data['type'])

    plt.show()

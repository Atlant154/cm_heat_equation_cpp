from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import math
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import ast

fig = plt.figure()
ax = fig.gca(projection='3d')

X_lin = np.linspace(0, 1, 100)
Y_lin = np.linspace(0, 1, 5)

X, Y = np.meshgrid(X_lin, Y_lin)

def exact_solution(X, Y):
    return - X ** 4 + Y * X + Y ** 2 - Y * np.exp(X)

Z = exact_solution(X, Y)

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

ax.set_xlabel('X')
ax.set_ylabel('Time')
ax.set_zlabel('Exact solution')

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
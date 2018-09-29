from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import ast

with open('result.txt', 'r') as file:
    for string in file:
        data = ast.literal_eval(string.strip())

fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.array(np.arange(np.float64(data[0][0]), np.float64(data[0][1]), np.float64(data[0][2])))
Y = np.array(np.arange(np.float64(data[1][0]), np.float64(data[1][1]), np.float64(data[1][2])))
Z = np.array(data[2])
X, Y = np.meshgrid(X, Y)


# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=1, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.show()

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
import ast

with open('error.txt', 'r') as file:
    for string in file:
        data = ast.literal_eval(string.strip())

fig = plt.figure()
ax = fig.gca(projection='3d')

X_LEFT_BOUND = np.float64(data[0][0])
X_RIGHT_BOUND = np.float64(data[0][1])
H_num = np.float64(data[0][2])
X_LIN = np.linspace(X_LEFT_BOUND, X_RIGHT_BOUND, H_num)

TIME_LEFT_BOUND = np.float64(data[1][0])
TIME_RIGHT_BOUND = np.float64(data[1][1])
Time_layers_num = np.float64(data[1][2])
Y_LIN = np.linspace(TIME_LEFT_BOUND, TIME_RIGHT_BOUND, Time_layers_num)

# Make data.
X = X_LIN
Y = Y_LIN
print(Y.shape)
print(X)
print(Y.shape)
Z = np.array(data[2])
print(Z.shape)
X, Y = np.meshgrid(X, Y)
print(X.shape)


# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=1, antialiased=False)

# Customize the z axis.
ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

ax.set_xlabel('X')
ax.set_ylabel('Time')
ax.set_zlabel('Temperature')

plt.show()

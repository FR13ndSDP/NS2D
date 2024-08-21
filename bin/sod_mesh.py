import numpy as np
import matplotlib.pyplot as plt

Lx = 1.0
Ly = 0.1

nx = 240
ny = 80

x = np.zeros((nx, ny))
y = np.zeros((nx, ny))

for i in range(nx):
    for j in range(ny):
        x[i, j] = i * Lx/(nx-1)
        y[i, j] = j * Ly/(ny-1)


# plt.scatter(x,y)
# plt.show()

with open("mesh2d.dat", "w") as f:
    f.write("variables=x,y\n")
    f.write("zone i=" + str(nx) + " j=" + str(ny) + "\n")
    for j in range(ny):
        for i in range(nx):
            f.write(str(x[i, j]) + " " + str(y[i, j]) + "\n")


import numpy as np
import matplotlib.pyplot as plt
Rin = 0.5
Rout = 1.5
Rstar = 0.5
alpha = 0.2
nx = 80
ny = 240

x = np.zeros((nx, ny))
y = np.zeros((nx, ny))

c1 = np.arcsin((Rin-Rstar)/alpha)
c2 = np.arcsinh((Rout-Rstar)/alpha)

for i in range(nx):
    for j in range(ny):
        theta = np.pi + np.pi/(ny-1) * j
        Ri = Rstar + alpha * np.sinh(c1*(1-i/(nx-1)) +c2*i/(nx-1))
        x[i, j] = Ri * np.sin(theta)
        y[i, j] = Ri * np.cos(theta)


# plt.scatter(x,y)
# plt.show()

with open("mesh2d.dat", "w") as f:
    f.write("variables=x,y\n")
    f.write("zone i=" + str(nx) + " j=" + str(ny) + "\n")
    for i in range(nx):
        for j in range(ny):
            f.write(str(x[i, j]) + " " + str(y[i, j]) + "\n")


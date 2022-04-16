import numpy as np
import matplotlib.pyplot as plt

def read_data():
    filename = 'final.dat'
    f   = open(filename)
    f.readline()
    f.readline()
    f.readline()
    x   = np.zeros((mx-1, my-1))
    y   = np.zeros((mx-1, my-1))
    rho = np.zeros((mx-1, my-1))
    u   = np.zeros((mx-1, my-1))
    v   = np.zeros((mx-1, my-1))
    p  = np.zeros((mx-1, my-1))
    T   = np.zeros((mx-1, my-1))
    i = 0; j = 0
    for line in f:
        x[i,j]   = eval(line.split()[0])
        y[i,j]   = eval(line.split()[1])
        rho[i,j] = eval(line.split()[2])
        u[i,j]   = eval(line.split()[3])
        v[i,j]   = eval(line.split()[4])
        p[i,j]  = eval(line.split()[5])
        T[i,j]   = eval(line.split()[6])

        if (i == mx-2):
            j += 1; i = 0
        else:
            i += 1
        if (j == my-1): break

    f.close()
    return (x,y,rho,u,v,p, T)
    
# grid points
mx = 80; my = 240 

gamma = 1.4
# read results
(x,y,rho,u,v,p, T) = read_data()

fig = plt.figure()
ax  = fig.add_subplot(111)

contour = ax.contourf(x,y,u, cmap=plt.cm.rainbow)
plt.colorbar(contour, orientation='vertical')
plt.axis("equal")
plt.show()
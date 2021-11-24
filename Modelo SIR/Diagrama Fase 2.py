# Código para graficar el comportamiento del diagrama de fase

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import pylab as pl
beta = 0.03
gamma = 0.1
N = 3.0
# The 3-dimensional nonlinear system.
def dx_dt(x, t):
 return [-beta*x[0]*x[1], beta*x[0]*x[1]-gamma*x[1], gamma*x[1]]
# Trajectories in forward time.
fig=plt.figure(figsize=(8, 6))
ts = np.linspace(0, 50000, 50000,dtype=int)
inc = np.linspace(0,2,15)
for s in inc:
 x0 = [N-s, s, 0]
 xs = odeint(dx_dt, x0, ts)
 plt.plot(xs[:,0], xs[:,1], "r-")
# plot I=N-S-R
x=np.linspace(0,3,25)
plt.plot(x,N-x,"k-")
# Label the axes and set fontsizes.
plt.xlabel("S", fontsize=15)
plt.ylabel("I", fontsize=15)
plt.tick_params(labelsize=15)
plt.xlim(0, 3)
plt.ylim(0, 3);
# Plot the vectorfield.
S, I = np.mgrid[0:3:20j, 0:3:20j]
u=-beta * S * I
v=beta * S * I - gamma * I
w=gamma * I
pl.quiver(S, I, u, v, color = 'b')
plt.savefig("Fase2.png")
plt.show()


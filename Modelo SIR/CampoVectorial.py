# Código que grafica el comportamiento de las soluciones del sistema

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline
from scipy.integrate import odeint
beta = 0.03
gamma = 0.1
N = 7.0
# The 3-dimensional nonlinear system.
def dx_dt(x, t):
 return [-beta*x[0]*x[1], beta*x[0]*x[1]-gamma*x[1], gamma*x[1]]
# Trajectories in forward time.
fig=plt.figure(figsize=(10, 8))
ax = fig.gca(projection='3d')
ts = np.linspace(0, 50000, 50000,dtype=int)
inc = np.linspace(0,3,15)
for s in inc:
 x0 = [N-s, s, 0]
 xs = odeint(dx_dt, x0, ts)
 ax.plot(xs[:,0], xs[:,1],xs[:,2], "r-")
# Plot the vectorfield.
S, I, R = np.mgrid[0:8:5j, 0:2:5j, 0:6:5j]
u=-beta * S * I
v=beta * S * I - gamma * I
w=gamma * I
ax.quiver(S,I,R,u,v,w,length=5, color = 'b')
ax.set_xlabel('S', fontsize=15)
ax.set_ylabel('I', fontsize=15)
ax.set_zlabel('R', fontsize=15)
plt.tick_params(labelsize=15)
ax.set_title('SIR Model', fontsize=15)
ax.xaxis.label.set_color('blue')
ax.yaxis.label.set_color('blue')
ax.zaxis.label.set_color('blue')
plt.savefig("SIR3.png")
plt.show()


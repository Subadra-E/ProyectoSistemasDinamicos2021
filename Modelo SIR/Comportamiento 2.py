# Código del comportamiento con un mejor ajuste

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# SIR model parameters and initial conditions
gamma, beta = 0.1, 0.03
S0, I0, R0 = 4.0-1.0, 1.0, 0.0
# Maximum time point and total number of time points
tmax, n = 50000, 100000
def SIR(X, t, gamma, beta):
 """SIR"""
 S, I, R = X
 dS = - beta * S * I
 dI = beta * S * I - gamma * I
 dR = gamma * I
 return dS, dI, dR
# Integrate the SIR model equations on the time grid t.
t = np.linspace(0, tmax, n)
f = odeint(SIR, (S0, I0, R0), t, args=(gamma, beta))
S, I, R = f.T
# Plot the SIR model using a Matplotlib 3D projection.
fig=plt.figure(figsize=(8, 6))
ax = Axes3D(fig)
ax.plot(S, I, R, 'r-', lw=2.0)
ax.set_xlabel('S', fontsize=15)
ax.set_ylabel('I', fontsize=15)
ax.set_zlabel('R', fontsize=15)
plt.tick_params(labelsize=15)
ax.set_title('SIR Model', fontsize=15)
ax.xaxis.label.set_color('blue')
ax.yaxis.label.set_color('blue')
ax.zaxis.label.set_color('blue')
plt.savefig("SIR1.png")
plt.show()


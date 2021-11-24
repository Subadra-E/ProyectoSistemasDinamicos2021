# Código para comparar curvas de infectados humanos


import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

N1 = 863689.0 #Population of humans
N2 = 1.0E05 #Population of mosquitoes
vh, vm, bhm, gh, bmh, uh, um  = 9.9982E-04, 8.8223E04, 8.8382E-05, 0.0991, 0.00362588, 0.01*7, 4.3E-05*7
Sh0, Ih0, Rh0, Sm0, Im0 = (N1-3.0)/N1, 3.0/N1, 0.0, N2-1, 1.0 
# Maximum time point and total number of time points
tmax, n = 53, 100000
def SIR1(X, t, vh, vm, bhm, gh, bmh, uh, um):
 """SIR with vectors"""
 Sh, Ih, Rh, Sm, Im = X
 dSh = vh - bhm * Im * Sh - uh * Sh
 dIh = bhm * Im * Sh - uh * Ih - gh * Ih
 dRh = -uh * Rh + gh * Ih
 dSm = vm - bmh * Ih * Sm - um * Sm
 dIm = bmh * Ih * Sm - um * Im
 return dSh, dIh, dRh, dSm, dIm
# Integrate the SIR model equations on the time grid t.
t = np.linspace(0, tmax, n)
f = odeint(SIR1, (Sh0, Ih0, Rh0, Sm0, Im0), t, args=(vh, vm, bhm, gh, bmh, uh, um))
Sh, Ih, Rh, Sm, Im = f.T
# Plot the SIR model using a Matplotlib 3D projection.
fig=plt.figure(figsize=(8, 6))
ax = Axes3D(fig)
ax.plot(Sh, Ih, Rh, 'b-', lw=2.0)
ax.set_xlabel('Sh', fontsize=15)
ax.set_ylabel('Ih', fontsize=15)
ax.set_zlabel('Rh', fontsize=15)
plt.tick_params(labelsize=15)
ax.set_title('SIR Model with vectors', fontsize=15)
ax.xaxis.label.set_color('red')
ax.yaxis.label.set_color('red')
ax.zaxis.label.set_color('red')
plt.savefig("SIRvectorsQtz.png")
plt.show()


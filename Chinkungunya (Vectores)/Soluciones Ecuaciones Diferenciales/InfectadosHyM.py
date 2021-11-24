# Código de comparación infectados humanos y mosquitos

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
# SIR model parameters and initial conditions
N1 = 3353951.0 #Population of humans
N2 = 1.0E05 #Population of mosquitoes
vh, vm, bhm, gh, bmh, uh, um  = 1.69E-06, 1.0E04, 1.62E-04, 0.08, 0.01, 0.01*7, 4.3E-05*7
Sh0, Ih0, Rh0, Sm0, Im0 = (N1-1.0)/N1, 1.0/N1, 0.0, N2-1, 1.0 
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
# Plot the SIR model
#SIR Humans
fig=plt.figure(figsize=(8, 6))
#plt.plot(t, Sh, 'g-', lw=2.0, label='Susceptible')
plt.plot(t, Ih, 'r-', lw=2.0, label='Infected')
plt.ylabel('Number of cases', fontsize=15)
plt.xlabel('Time [weeks]', fontsize=15)
plt.tick_params(labelsize=15)
plt.title('SIR Model with vectors', fontsize=15)
lgd = plt.legend(bbox_to_anchor=(1.01,0.65), loc="center left",fontsize=15)
plt.tight_layout()
plt.savefig("SIRvectors.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
#SIR Mosquitoes
fig1=plt.figure(figsize=(8, 6))
#plt.plot(t, Sm/N2, 'g--', lw=2.0, label='Susceptible mosquitoes')
plt.plot(t, Im/N2, 'r--', lw=2.0, label='Infected mosquitoes')
plt.ylabel('Number of cases', fontsize=15)
plt.xlabel('Time [weeks]', fontsize=15)
plt.tick_params(labelsize=15)
plt.title('SIR Model with vectors', fontsize=15)
lgd1 = plt.legend(bbox_to_anchor=(1.01,0.65), loc="center left",fontsize=15)
plt.tight_layout()
plt.savefig("SIRvectors1.png", bbox_extra_artists=(lgd1,), bbox_inches='tight')
plt.show()
#Infected humans and mosquitoes
fig2=plt.figure(figsize=(8, 6))
#plt.plot(t, Sh, 'g-', lw=2.0, label='Susceptible')
plt.plot(t, Ih, 'r-', lw=2.0, label='Infected')
plt.plot(t, Im/N2, 'r--', lw=2.0, label='Infected mosquitoes')
plt.ylabel('Number of cases', fontsize=15)
plt.xlabel('Time [weeks]', fontsize=15)
plt.tick_params(labelsize=15)
plt.title('SIR Model with vectors', fontsize=15)
lgd2 = plt.legend(bbox_to_anchor=(1.01,0.65), loc="center left",fontsize=15)
plt.tight_layout()
plt.savefig("Infecteshm.png", bbox_extra_artists=(lgd2,), bbox_inches='tight')


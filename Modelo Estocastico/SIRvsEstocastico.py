# Comparación de modelos


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint


def f(y, t, paras):
    """
    Your system of differential equations
    """

    x1 = y[0] #S
    x2 = y[1] #I
    x3 = y[2] #R

    try:
        k0 = paras['k0'].value #beta (transmission rate)
        k1 = paras['k1'].value #gamma (recovery rate)

    except KeyError:
        k0, k1 = paras
    # the model equations
    f0 = -k0 * x1 * x2 #dS
    f1 = k0 * x1 * x2- k1 * x2 #dI
    f2 = k1 * x2 #dR
    return [f0, f1, f2]


def g(t, x0, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(f, x0, t, args=(paras,))
    return x


def residual(paras, t, data):

    """
    compute the residual between actual data and fitted data
    """

    x0 = paras['x10'].value, paras['x20'].value, paras['x30'].value
    model = g(t, x0, paras)

    # you only have data for one of your variables
    x2_model = model[:, 1]
    return (x2_model - data).ravel()


# initial conditions
N = 10000. #Total population
x10 = 9995/1579.0 # S(0)
x20 = 5/1579.0  # I(0)
x30 = 0.  # R(0)
y0 = [x10, x20, x30]

#Here the stochastic simulation data is substituted 

# measured data
t_measured = list(np.loadtxt("MedianStats.txt")[:, 0]) 
data = list(np.loadtxt("MedianStats.txt")[:, 2])
median_measured = list(map(lambda x: x / 1579.0, data))
#lower_measured = np.loadtxt("LowerStats.txt")[:, 2]
#higher_measured = np.loadtxt("HigherStats.txt")[:, 2]

plt.figure()
fig=plt.figure(figsize=(8, 6))
plt.scatter(t_measured, median_measured, marker='o', color='r', label='measured data', s=30)

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('x10', value=x10, vary=False)
params.add('x20', value=x20, vary=False)
params.add('x30', value=x30, vary=False)
params.add('k0', value=0.03, min=0.0001, max=2.)
params.add('k1', value=0.1, min=0.0001, max=2.)

# fit model
result = minimize(residual, params, args=(t_measured, median_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(np.linspace(0., 149., 100), y0, result.params)

# plot fitted data
plt.plot(np.linspace(0., 149., 100), data_fitted[:, 1], '-', linewidth=2.5, color='blue', label='fitted data')
plt.legend()
plt.xlim([0, max(t_measured)])
plt.ylim([0, 1.1 * max(data_fitted[:, 1])])
plt.xlabel("Time steps [days]",size=15)
plt.ylabel("Number of cases",size=15)
lgd = plt.legend(bbox_to_anchor=(1.01,0.65), loc="center left",fontsize=15)
plt.tight_layout()
plt.savefig("comparacion2v1.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
# display fitted statistics
report_fit(result)

plt.show()




import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint
N1  = 863689. #Population of Quetzaltenango
N2 = 1.0E06 #Population of mosquitoes

def f(y, t, paras):
    """
    Your system of differential equations
    """
    x1 = y[0]                   #Sh
    x2 = y[1]                   #Ih
    x3 = y[2]                   #Rh
    x4 = y[3]                   #Sm
    x5 = y[4]                   #Im

    try:
        vh = paras['vh'].value
        vm = paras['vm'].value
        bhm = paras['bhm'].value
        gh = paras['gh'].value
        bmh = paras['bmh'].value
        uh = paras['uh'].value
        um = paras['um'].value
    except KeyError:
      vh, vm, bhm, gh, bmh, uh, um  = paras
    # the model equations
    f1 = vh - bhm * x1 * x5 - uh * x1 #dSh
    f2 = bhm * x1 * x5 - uh * x2 - gh * x2 #dIh
    f3 = -uh * x3 + gh * x2 #dRh
    f4 = vm - bmh * x2 * x4 - um * x4 #dSm
    f5 = bmh * x2 * x4 - um * x5 #dIm

    return [f1, f2, f3, f4, f5]

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

    x0 = paras['x10'].value, paras['x20'].value, paras['x30'].value, paras['x40'].value, paras['x50'].value
    model = g(t, x0, paras)

    # you only have data for one of your variables
    x2_model = model[:, 1]
    return (x2_model - data).ravel()


# initial conditions
x10 = (N1-3.0)/N1 #Sh(0)
x20 = 3.0/N1 #Ih(0)
x30 = 0.0/N1 #Rh(0)
x40 = N2-1.0 #Sm(0)
x50 = 1.0 #Im(0)
y0 = [x10, x20, x30, x40, x50]

# measured data
data = list(np.loadtxt("Quetzaltenango.txt"))
x2_measured = list(map(lambda x: x / 181.0, data))

t_measured = np.linspace(1, 52,len(x2_measured))
plt.figure(figsize=(8, 6))
plt.scatter(t_measured, x2_measured, marker='o', color='b', label='measured data', s=20)

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('x10', value=x10, vary=False)
params.add('x20', value=x20, vary=False)
params.add('x30', value=x30, vary=False)
params.add('x40', value=x40, vary=False)
params.add('x50', value=x50, vary=False)
params.add('vh', value=1.69E-06, min=1.0E-07, max=1.0E-03)#lo moví, muy sensible
params.add('vm', value=1.0E04, min=1.0E03, max=1.0E06)
params.add('vh', value=1.69E-06, min=1.0E-07, max=1.0E-03)
params.add('bhm', value=1.62E-04, min=1.0E-06, max=1.0E-02 )
params.add('gh', value=0.08, min=1E-03, max=1)
params.add('bmh', value=2.18E-02, min=1.0E-05, max=1.0E-02)
params.add('uh', value=4.3E-05*7 , vary=False)
params.add('um', value=0.01*7, vary=False)

# fit model
result = minimize(residual, params, args=(t_measured, x2_measured), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(t_measured, y0, result.params)

# plot fitted data
plt.plot(t_measured, data_fitted[:, 1], '-', linewidth=2, color='red', label='fitted data')
plt.legend()
plt.xlim([0, max(t_measured)])
plt.ylim([0, 1.1 * max(data_fitted[:, 1])])
plt.xlabel("Time steps [weeks]",size=15)
plt.ylabel("Number of cases",size=15)
lgd = plt.legend(bbox_to_anchor=(1.01,0.65), loc="center left",fontsize=15)
plt.tight_layout()
plt.savefig("QuetzaltenangoChinkungunya.png", bbox_extra_artists=(lgd,), bbox_inches='tight')
# display fitted statistics
report_fit(result)
plt.show()


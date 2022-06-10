import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('font', size = 16)                 # Use big fonts...
plt.rcParams['figure.figsize'] = (10, 8)         # ... and big plots

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm


import CoolProp.CoolProp as cp
import numpy as np

flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)

c_vap = cp.PropsSI('H', 'P', 4.7*101325, 'Q', 1, 'Xenon') - cp.PropsSI('H', 'P', 4.7*101325, 'Q', 0, 'Xenon')
print(c_vap, 'J/kg')

dqdl = 6.4  # W/m
mdot = 2 * 2.9 / 60  # kg/s
print(c_vap / (dqdl/mdot), 'm')
print(c_vap / ((dqdl/mdot) + 9.8), 'm')

def get_h(mdot_lpm, pbar=4., Tc=-90.0, reduce=True):
    venturi_coeff = 92.16  # (lpm^2)/bar
    h = 0.
    dh = 0.01  # height step
    pa_per_bar = 101325
    dqdl = 6.4  # W/m
    mdot = mdot_lpm * 2.9 / 60  # kg/s
    venturi_reduction = (mdot_lpm**2.0)/venturi_coeff
    fcv286_cv = 2.35
    fcv286_reduction = (mdot_lpm/fcv286_cv)**2.0/72.55
    #print(venturi_reduction, venturi_reduction)
    
    pbar = pbar - venturi_reduction - fcv286_reduction  # 
    p = pbar * pa_per_bar
    T = Tc + 273.15
    rho = cp.PropsSI('D', 'P', p, 'T', T, 'Xenon')
    enth = cp.PropsSI('H', 'P', p, 'T', T, 'Xenon')
    while 1:
        h += dh
        if h>100:
            break
        p -= 9.8*dh*rho  # pressure decreases due to height
        enth += (dqdl * dh / mdot) - (9.8 * dh)  # enth increases due to heat input, decreases due to height drop
        T = cp.PropsSI('T', 'H', enth, 'P', p, 'Xenon')  # temperature at new enthalpy and pressure
        Tsat = cp.PropsSI('T', 'P', p, 'Q', 0, 'Xenon')  # Tsat at same pressure
        #print(p, enth, Tsat - T)
        if (Tsat-T) < 0.1:
            break
    return h

def get_l(mdot_lpm, pbar=4., Tc=-90.0):
    h = 0.
    dh = 0.1  # height step
    pa_per_bar = 101325
    dqdl = 6.4  # W/m
    mdot = mdot_lpm * 2.9 / 60  # kg/s
    p = pbar*pa_per_bar
    T = Tc+273.15
    rho = cp.PropsSI('D', 'P', p, 'T', T, 'Xenon')
    enth = cp.PropsSI('H', 'P', p, 'T', T, 'Xenon')
    while 1:
        h += dh
        if h>100:
            break
        enth += (dqdl * dh / mdot)  # enth increases due to heat input, decreases due to height drop
        T = cp.PropsSI('T', 'H', enth, 'P', p, 'Xenon')  # temperature at new enthalpy and pressure
        Tsat = cp.PropsSI('T', 'P', p, 'Q', 0, 'Xenon')  # Tsat at same pressure
        #print(p, enth, Tsat - T)
        if (Tsat-T) < 0.05:
            break
    return h

minflow = flow_func(61.5)
minheight = get_h(minflow, pbar=4.05, Tc=-90.2)

flows = {}
hs = {}
for T in [-90., -95.]:
    for pbar in [4.0, 4.5, 5.0]:
        flows[(T, pbar)] = np.linspace(0.1, 3.0, 100)
        hs[(T, pbar)] = [get_h(md, pbar=pbar, Tc=T) for md in flows[(T, pbar)]]

for T, pbar in (
    (-95., 5.0),
    (-95., 4.5),
    (-90., 5.0),
    (-95., 4.0),
    (-90., 4.5),
    (-90., 4.0),
):
    plt.plot(flows[(T, pbar)], hs[(T, pbar)], label='%.1f bar, %d C' % (pbar, int(T)))
plt.plot([minflow], [minheight], 'kx', markersize=10)
plt.axhline(7.25, color='k', linestyle='--')
plt.xlabel('LXe Flow [liter/min]')
plt.ylabel('Saturation Height [m]')
plt.annotate('4.05 bar, -90.2 C', (1.65, 5.25))
plt.annotate('Tubing Max Height', (2.1, 7.5))
plt.legend(loc=4)
plt.grid(alpha=0.3)
plt.xlim([0, 3])
if True:
    plt.savefig('plots/satheight.png', dpi=200)
    plt.savefig('plots/satheight.pdf')
plt.show()

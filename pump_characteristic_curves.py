## Code taken from Guillaume

import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
matplotlib.rc('font', size = 16)                 # Use big fonts...
plt.rcParams['figure.figsize'] = (10, 8)         # ... and big plots

import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

import time, calendar

save = True

## P-271

csv_files = [
    'data/xenonnt-lxe_purification-pump_characterization-run00_commissioning-20200131.csv',
    'data/xenonnt-lxe_purification-pump_characterization-run00_commissioning-20200201.csv'
]

time_intervals = {}

# 38 Hz
time_intervals[38] = [
        ['01/31/20 09:00:00', '01/31/20 09:20:00'],
        ['01/31/20 09:42:00', '01/31/20 09:52:00'],
        ['01/31/20 10:01:00', '01/31/20 10:11:00'],
        ['01/31/20 10:17:00', '01/31/20 10:27:00'],
        ['01/31/20 14:26:00', '01/31/20 14:36:00'],
        ['01/31/20 14:39:00', '01/31/20 14:49:00'],
        ['01/31/20 14:52:30', '01/31/20 15:02:30'],
        ['01/31/20 15:08:00', '01/31/20 15:18:00'],
        ['01/31/20 15:23:00', '01/31/20 15:33:00'], # maybe low flow
        ['01/31/20 15:40:00', '01/31/20 15:50:00'], # maybe low flow
        ['01/31/20 16:04:00', '01/31/20 16:14:00'], # maybe low flow
        ['01/31/20 16:23:00', '01/31/20 16:33:00'], # maybe low flow
        ['01/31/20 16:50:00', '01/31/20 17:00:00'],
        ['02/01/20 11:11:00', '02/01/20 11:25:00'],
        ['02/01/20 11:43:00', '02/01/20 11:58:00'],
        ['02/01/20 13:00:00', '02/01/20 13:30:00'],
        ['02/01/20 16:30:00', '02/01/20 16:45:00']
        ]
# 40 Hz
time_intervals[40] = [
        ['01/31/20 10:37:00', '01/31/20 10:47:00'],
        ['01/31/20 10:57:00', '01/31/20 11:07:00'],
        ['01/31/20 11:11:00', '01/31/20 11:21:00'],
        ['01/31/20 11:27:00', '01/31/20 11:37:00'],
        ['01/31/20 11:43:00', '01/31/20 11:53:00'],
        ['02/01/20 17:40:00', '02/01/20 18:00:00']
        ]
# 42 Hz
time_intervals[42] = [
        ['01/31/20 12:02:00', '01/31/20 12:12:00'],
        ['01/31/20 12:18:00', '01/31/20 12:28:00'],
        ['01/31/20 12:33:00', '01/31/20 12:43:00'],
        ['01/31/20 12:48:00', '01/31/20 12:58:00'],
        ['02/01/20 14:30:00', '02/01/20 14:45:00']
        ]
# 44 Hz
time_intervals[44] = [
        ['01/31/20 13:15:00', '01/31/20 13:20:00'],
        ['01/31/20 13:24:00', '01/31/20 13:34:00'],
        ['01/31/20 13:38:00', '01/31/20 13:48:00'],
        ['01/31/20 13:53:00', '01/31/20 14:03:00'],
        ['01/31/20 14:07:00', '01/31/20 14:17:00'],
        ['02/01/20 15:00:00', '02/01/20 15:20:00']
        ]

pressure_offsets = {'pt101':0.0000, 'pt271':+0.0380, 'pt273':+0.0296, 'pt275':+0.0036, 'pt286':+0.0237, 'pt287':+0.0511}

hv277_cv = 4.2
hv278_cv = 4.2
fv272_cv = 4.2
fv274_cv = 4.2
pcv287_cv = 2.35

to_unixtime = lambda t: calendar.timegm(time.strptime(t[:17], '%m/%d/%y %H:%M:%S'))-1*3600
flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)

dataframes_list = [pd.read_csv(f, converters={'LNGS Time': lambda t: calendar.timegm(time.strptime(t, '%Y-%m-%d %H:%M:%S'))-1*3600}) for f in csv_files]
sensors_dataframe = pd.concat(dataframes_list)

sensors_dataframe = sensors_dataframe.astype(np.float64)
unixtime_vals_array = sensors_dataframe['LNGS Time'].values.astype(np.float64)
pt271_vals_array = sensors_dataframe['PUR_PT271_PHEOUTLET.PI'].values.astype(np.float64)
pt275_vals_array = sensors_dataframe['PUR_PT275_PPUMPOUTLET.PI'].values.astype(np.float64)

dpt288_vals_array = sensors_dataframe['PUR_DPT288_PMON_A.PI'].values.astype(np.float64)

fcv286_vals_array = sensors_dataframe['PUR_FCV286SET_AO.PI'].values.astype(np.float64)
pcv287_vals_array = sensors_dataframe['PUR_PCV287SET_AO.PI'].values.astype(np.float64)

p271_fmon_vals_array = sensors_dataframe['PUR_P271_P271_PUB_FMON.PI'].values

pump_frequencies = time_intervals.keys()

## aggregate
flows = {}
deltaps = {}
for pump_frequency in pump_frequencies:
    flow_vals = []
    deltap_vals = []
    for time_interval in time_intervals[pump_frequency]:
        t0, t1 = tuple(map(to_unixtime, time_interval))
        interval_array = (unixtime_vals_array >= t0) & (unixtime_vals_array <= t1)

        fcv286_opening = fcv286_vals_array[interval_array][0]
        if not np.all(fcv286_vals_array[interval_array] == fcv286_opening):
                print('warning: FCV286 opening not constant over entire interval')

        pcv287_opening = pcv287_vals_array[interval_array][0]
        if not np.all(pcv287_vals_array[interval_array] == pcv287_opening):
                print('warning: PCV287 opening not constant over entire interval')

        pt271_val = np.mean(pt271_vals_array[interval_array])-pressure_offsets['pt271']
        pt275_val = np.mean(pt275_vals_array[interval_array])-pressure_offsets['pt275']

        dpt288_val = np.mean(dpt288_vals_array[interval_array])
        flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)
        dpt288_flow_val = flow_func(dpt288_val)

        p271_fmon_val = p271_fmon_vals_array[interval_array][0]
        if not np.all(p271_fmon_vals_array[interval_array] == p271_fmon_val):
                print('warning: pump frequency not constant over entire interval')

        C1 = 72.55*(pt275_val-pt271_val)
        C2 = 1./fv272_cv**2
        C3 = 1./(0.01*pcv287_opening*pcv287_cv)**2

        discr = np.sqrt((2.*C3*dpt288_flow_val)**2 - 4.*(C2+C3)*(C3*dpt288_flow_val**2-C1))

        flow_a = (2.*C3*dpt288_flow_val + discr)/(2.*(C2+C3))
        flow_b = (2.*C3*dpt288_flow_val - discr)/(2.*(C2+C3))

        flow = flow_a

        deltap = ((1./hv277_cv**2 + 1./hv277_cv**2)*flow**2 + (flow-dpt288_flow_val)**2/(0.01*pcv287_opening*pcv287_cv)**2)/72.55

        flow_vals.append(flow)
        deltap_vals.append(deltap)
    flows[pump_frequency] = flow_vals
    deltaps[pump_frequency] = deltap_vals


## model
from scipy.optimize import curve_fit, minimize

def deltap_for_curvefit(
        x0,
        a0,
        a1,
        b0,
        b1,
        b2,
        c0
    ):
    freq = 44.0
    mdot = x0
    deltap = ((a0 + (a1*freq))*(b0 + (b1*mdot) + (b2*(mdot**2.0)))) + c0
    return deltap

def deltap_model_fixed_p(
        mdot,
        freq,
        a0,
        a1,
        b0,
        b1,
        b2,
        c0
    ):
    deltap = ((a0 + (a1*freq))*(b0 + (b1*mdot) + (b2*(mdot**2.0)))) + c0
    return deltap

def get_chi2(data, model):
    # both array-like
    chi2 = np.sum((data - model)**2.0/(2*model**2.0))
    return chi2

def get_neg_lnl(data, model, uncertainty):
    # both arrays equal size
    chi2 = np.sum((data - model)**2.0/(2*uncertainty**2.0))
    normalization = np.sum(-1*np.log(np.sqrt(2*np.pi)) - np.log(uncertainty))
    lnl = normalization - chi2
    return -1.0*lnl

def full_lnl(pars, data):
    neg_lnl = 0
    for row in data:
        model = deltap_model_fixed_p(
                row[1],
                row[0],
                *pars
            )
        neg_lnl += get_chi2(
            row[2],
            model
        )
    return neg_lnl

# guess
a0, a1 = (0.2,   1.9)
b0, b1, b2 = (3.426872e-2,   1.83331873e-03, -1.53541442e-04)
c0 = -1.2

freqs = np.array(list(flows.keys()))
print(freqs)

data = [
    (
        freq,
        np.array(flows[freq]),
        np.array(deltaps[freq])
    )
    for freq in freqs
]

res = minimize(
    full_lnl,
    (a0, a1, b0, b1, b2, c0),
    args=(data),
    options={'maxiter': 100000000},
    method='TNC',
    bounds=(
        (None, None),
        (None, None),
        (None, None),
        (None, None),
        (-np.inf, 0),
        (None, None),
    )
)
print(res.x)

xs = np.linspace(3, 8, 100)
guess = False
colors = (np.array(freqs) - freqs[0])/(freqs[-1] - freqs[0])
cmap = cm.get_cmap('viridis', len(freqs))
textpos = [
    (4.95, 1.43),
    (4.65, 1.63),
    (4.85, 1.79),
    (5, 1.96),
]
for i, freq in enumerate(freqs):
    plt.plot(
        flows[freq], 
        deltaps[freq],
        marker='.',
        color=cmap(colors[i]),
        lw=0,
        label='%i Hz' % int(freq)
    )
    pars = res.x
    if guess:
        pars = (a0, a1, b0, b1, b2, c0)
    ys = deltap_model_fixed_p(
        xs,
        freq,
        *pars
    )
    plt.plot(xs, ys, linestyle='--', color=cmap(colors[i]))
    plt.annotate('%i Hz' % int(freq), textpos[i], color=cmap(colors[i]))
#plt.legend()
plt.xlabel('Flow [l/min]')
plt.ylabel('Pump dP [bar]')
if save:
    name = 'plots/P271_characterization'
    plt.savefig(name+'.png', dpi=300)
    plt.savefig(name+'.pdf')
names = ['a0', 'a1', 'b0', 'b1', 'b2', 'c0']
for i, v in enumerate(res.x):
    print('|', names[i], '|', '%.2e' % v, '|')
plt.show()


### P-272
csv_files = [
    'data/data-from -02_06_2020-08 00 00-to -02_06_2020-20 00 00.csv',
    'data/data-from -02_08_2020-05 00 00-to -02_08_2020-22 00 00.csv',
    'data/data-from -12_06_2021-08 00 00-to -12_31_2021-20 00 00.csv',
]
time_intervals = {}

# 31 Hz
time_intervals[31] = [
    ['12/16/21 23:00:00', '12/19/21 16:00:00'],
    ['12/21/21 22:00:00', '12/28/21 21:00:00'],
]
# 32 Hz
time_intervals[32] = [
    ['12/06/21 21:00:00', '12/10/21 21:00:00'],
    ['12/20/21 19:00:00', '12/21/21 09:00:00'],
]
# 33 Hz
time_intervals[33] = [
    ['02/06/20 12:00:00', '02/06/20 12:10:00'],
    ['02/06/20 17:05:00', '02/06/20 18:29:00'],
    ['02/08/20 06:20:00', '02/08/20 08:31:00'],
    ['02/08/20 10:37:00', '02/08/20 10:54:00'],
    ['02/08/20 10:59:00', '02/08/20 11:45:00'],
    ['02/08/20 14:55:00', '02/08/20 15:20:00'],
]
# 34 Hz
time_intervals[34] = [
    ['02/06/20 12:24:00', '02/06/20 12:32:00'],
    ['02/06/20 12:36:00', '02/06/20 13:48:00'],
    ['02/06/20 15:27:00', '02/06/20 16:50:00'],
    ['02/06/20 18:50:00', '02/06/20 19:47:00'],
    ['02/08/20 14:07:00', '02/08/20 14:43:00'],
    ['02/08/20 16:30:00', '02/08/20 16:55:00'],
    ['02/08/20 17:06:00', '02/08/20 20:00:00'],
    ['12/11/21 18:00:00', '12/15/21 04:00:00'],
]
# 35 Hz
time_intervals[35] = [
    ['02/08/20 12:02:00', '02/08/20 12:42:00'],
]
# 36 Hz
time_intervals[36] = [
    ['02/06/20 14:02:00', '02/06/20 14:43:00'],
    ['02/06/20 14:53:00', '02/06/20 15:20:00'],
]
# 38 Hz
time_intervals[38] = [
    ['02/08/20 13:04:00', '02/08/20 13:43:00'],
]

dataframes_list = [pd.read_csv(f, converters={'LNGS Time': lambda t: calendar.timegm(time.strptime(t, '%Y-%m-%d %H:%M:%S'))-1*3600}) for f in csv_files]
sensors_dataframe = pd.concat(dataframes_list)

unixtime_vals_array = sensors_dataframe['LNGS Time'].values.astype(np.float64)
pt271_vals_array = sensors_dataframe['PUR_PT271_PHEOUTLET.PI'].values.astype(np.float64)
pt276_vals_array = sensors_dataframe['PUR_PT276_PPUMPOUTLET.PI'].values.astype(np.float64)

dpt288_vals_array = sensors_dataframe['PUR_DPT288_PMON_A.PI'].values.astype(np.float64)

fcv286_vals_array = sensors_dataframe['PUR_FCV286SET_AO.PI'].values.astype(np.float64)
pcv287_vals_array = sensors_dataframe['PUR_PCV287SET_AO.PI'].values.astype(np.float64)

p271_fmon_vals_array = sensors_dataframe['PUR_P271_P271_PUB_FMON.PI'].values

pump_frequencies = time_intervals.keys()

flows = {}
deltaps = {}
print('flow, dpt288_flow_val, fcv286_opening, pcv287_opening')
for pump_frequency in pump_frequencies:
    #print(pump_frequency)
    flow_vals = []
    deltap_vals = []
    if pump_frequency==33:
        fig = plt.figure()
    for time_interval in time_intervals[pump_frequency]:
        if pump_frequency==33:
            print(time_interval)
        t0, t1 = tuple(map(to_unixtime, time_interval))
        interval_array = (unixtime_vals_array >= t0) & (unixtime_vals_array <= t1)

        fcv286_opening = fcv286_vals_array[interval_array][0]
        if not np.all(fcv286_vals_array[interval_array] == fcv286_opening):
                print('warning: FCV286 opening not constant over entire interval')

        pcv287_opening = pcv287_vals_array[interval_array][0]
        if not np.all(pcv287_vals_array[interval_array] == pcv287_opening):
                print('warning: PCV287 opening not constant over entire interval')

        pt271_val = np.mean(pt271_vals_array[interval_array])-pressure_offsets['pt271']
        pt276_val = np.mean(pt276_vals_array[interval_array])-pressure_offsets.get('pt286', 0.0)

        dpt288_val = np.mean(dpt288_vals_array[interval_array])
        flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)
        dpt288_flow_val = flow_func(dpt288_val)

        p271_fmon_val = p271_fmon_vals_array[interval_array][0]
        if not np.all(p271_fmon_vals_array[interval_array] == p271_fmon_val):
                print('warning: pump frequency not constant over entire interval')

        C1 = 72.55*(pt276_val-pt271_val)
        C2 = 1./fv272_cv**2
        C3 = 1./(0.01*pcv287_opening*pcv287_cv)**2

        discr = np.sqrt((2.*C3*dpt288_flow_val)**2 - 4.*(C2+C3)*(C3*dpt288_flow_val**2-C1))

        flow_a = (2.*C3*dpt288_flow_val + discr)/(2.*(C2+C3))
        flow_b = (2.*C3*dpt288_flow_val - discr)/(2.*(C2+C3))

        flow = flow_a

        deltap = ((1./hv277_cv**2 + 1./hv277_cv**2)*flow**2 + (flow-dpt288_flow_val)**2/(0.01*pcv287_opening*pcv287_cv)**2)/72.55
        flow_vals.append(flow)
        deltap_vals.append(deltap)
    flows[pump_frequency] = flow_vals
    deltaps[pump_frequency] = deltap_vals

# guess
a0, a1, b0, b1, b2, c0 = (
    0.3960738681121542, 
    0.6281057826493459, 
    0.16349687048634157, 
    0.000596546113237772, 
    -3.430201196170388e-05, 
    -1.7311328538521062,
)
freqs = np.array(list(flows.keys()))
print(freqs)

data = [
    (
        freq,
        np.array(flows[freq]),
        np.array(deltaps[freq])
    )
    for freq in freqs
]

res = minimize(
    full_lnl,
    (a0, a1, b0, b1, b2, c0),
    args=(data),
    options={'maxiter': 100000000},
    method='TNC',
    bounds=(
        (None, None),
        (None, None),
        (None, None),
        (None, None),
        (-np.inf, 0),
        (None, None),
    )
)
xs = np.linspace(3, 13, 100)
guess = False
colors = (np.array(freqs) - freqs[0])/(freqs[-1] - freqs[0])
cmap = cm.get_cmap('viridis', len(freqs))
textpos = [
    (9.6, 1.59),
    (9.7, 1.69),
    (9.8, 1.79),
    (9.9, 1.9),
    (10, 2.0),
    (10.1, 2.1),
    (10.3, 2.25),
]
for i, freq in enumerate(freqs):
    print(freq)
    print(flows[freq])
    plt.plot(
        flows[freq], 
        deltaps[freq],
        marker='.',
        color=cmap(colors[i]),
        lw=0,
        label='%i Hz' % int(freq)
    )
    pars = res.x
    if guess:
        pars = (a0, a1, b0, b1, b2, c0)
    ys = deltap_model_fixed_p(
        xs,
        freq,
        *pars
    )
    plt.plot(xs, ys, linestyle='--', color=cmap(colors[i]))
    plt.annotate('%i Hz' % int(freq), textpos[i], color=cmap(colors[i]))
plt.xlabel('Pump Flow [liter/min]')
plt.ylabel(r'Pump $\Delta P$ [bar]')
if save:
    name = 'plots/P272_characterization'
    plt.savefig(name+'.png', dpi=300)
    plt.savefig(name+'.pdf')
plt.grid('on')
names = ['a0', 'a1', 'b0', 'b1', 'b2', 'c0']
for i, v in enumerate(res.x):
    print('|', names[i], '|', '%.2e' % v, '|')
plt.show()

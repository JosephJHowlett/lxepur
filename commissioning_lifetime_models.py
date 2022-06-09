import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('font', size = 18)                 # Use big fonts...
plt.rcParams['figure.figsize'] = (10, 8)         # ... and big plots
plt.rcParams['timezone'] = 'Europe/Paris'

import pandas as pd
import datetime
import numpy as np
import matplotlib.dates as dates
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines

import datetime as dt
from scipy.signal import medfilt

def string_to_lngs_epoch(string, datetime_format='%m-%d-%y %H:%M:%S'):
    # Just putting things in LNGS time and making a unified date format
    # to make things simple.
    epoch = datetime.datetime(1970, 1, 1)
    return (datetime.datetime.strptime(string, datetime_format) - epoch).total_seconds()

def get_errs(elifetime, anodeheight, cathodeheight, drifttime):
    drifttime_uncertainty = 0.7 # us
    anodeheight_uncertainty = 0.2 # mV
    cathodeheight_uncertainty = 0.2 # mV
    return (
        elifetime*np.sqrt(
            (drifttime_uncertainty/drifttime)**2 
             + (anodeheight_uncertainty/(anodeheight*np.log(anodeheight/cathodeheight)))**2 
            + (cathodeheight_uncertainty/(cathodeheight*np.log(anodeheight/cathodeheight)))**2)
    )

flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)

use_scada = False

if use_scada:
    from straxen import SCADAInterface
    sc = SCADAInterface()
    start_time = '10-13-20 09:02:00'
    stop_time = '12-01-20 04:02:00'  # go till 12/1 and see

    var_names = {
        'lifetime': 'XE1T.PUR_LXEELIFETIME.PI',
        'anode_height': 'XE1T.PUR_LXEANODEPULSEHIGHT.PI',
        'cathode_height': 'XE1T.PUR_LXECATODEPULSEHIGHT.PI',
        'drift_time': 'XE1T.PUR_LXEEDRIFTTIME.PI',
    }

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
else:
    df = pd.read_csv('data/data_commissioning_lifetime_q5.csv', index_col='time UTC', parse_dates=True)

# clean
df_clean = df[
    ((df.lifetime < 9e3) | (df.index > dt.datetime(2020, 11, 8, tzinfo=dt.timezone.utc)))
    & (df.lifetime > 0)
]
df_clean = df_clean[df_clean.index > dt.datetime(2020, 10, 23, tzinfo=dt.timezone.utc)]
tdiff_days = (df_clean.index - df_clean.index[0]).total_seconds()/60./60./24.
err = get_errs(df_clean.lifetime, df_clean.anode_height, df_clean.cathode_height, df_clean.drift_time)
lt = medfilt(df_clean.lifetime, kernel_size=501)

## model
from scipy.integrate import odeint

def get_RHS(x, t, p, t1=np.inf):
    rhs = p[1]
    rhs += -1.0*x/p[0]
    if t>t1:
        rhs += -1.0*x/p[2]
    return rhs

conv_outg = lambda outg: outg * 257.0 * 2.9 / 1e9 * 32 / 131 * 1000.0  # conv-to-mass g/d
inv_outg = lambda outg: outg / 257.0 / 2.9 * 1e9 / 32 * 131 / 1000.0  # conv-to-us us^-1*l/d

density = 18.694  # kg/m^3 
density = 5.887  # g/SL
lxe_density = 2.9  # kg/l

lxe_mass = 8600.  # kg
outgassing = inv_outg(0.06e-3)  # g/d
outgassing = inv_outg(0.11e-3)  # g/d
flow_speed = 60.0
lxe_flow_speed = 2.0  # lpm

tmax = 40.0
x0 = 0.0142
t1 = 14.2

times = np.linspace(0, tmax, 10000)

tau = lxe_mass/density*1000/flow_speed*60/3600/24  # d
tau2 = lxe_mass/lxe_density/lxe_flow_speed*60/3600/24  # d
print('outgassing', 'tau', 'tau2')
print(outgassing, tau, tau2)
plt.figure(figsize=(12, 8))

plt.axvline(t1, color='C1', linestyle=':')
plt.annotate('Start\nLXe Purification', (14.7, 11e3), color='C1', ha='left')
plt.plot(tdiff_days, lt, 'k.', rasterized=True)
plt.fill_between(tdiff_days, lt-err, lt+err, color='slategrey', alpha=0.3, linewidth=0.2)

p = [tau, outgassing / (lxe_mass/lxe_density), tau2]
sol = odeint(
    get_RHS,
    x0,
    times,
    args=(p, np.inf)
)
plt.plot(times, 1.0/sol.T[0], 'C0--', lw=3)
p = [tau, outgassing / (lxe_mass/lxe_density), tau2]
sol = odeint(
    get_RHS,
    x0,
    times,
    args=(p, t1)
)
plt.plot(times, 1.0/sol.T[0], 'C0-', lw=3)
plt.yscale('log')
plt.grid(alpha=0.3)
plt.xlim([0, 40])
#plt.ylim([0, 6000])
plt.xlabel('Time [d]')
plt.ylabel('Electron Lifetime [$\mu$s]')
if True:
    plt.savefig('plots/q5_trend.pdf')
    plt.savefig('plots/q5_trend.png', dpi=200)
plt.show()



### ST707

if use_scada:
    start_time = '01-19-21 00:00:00'
    stop_time = '02-03-21 00:00:00'
    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
else:
    df = pd.read_csv('data/st707commissioning_full.csv', index_col='time UTC', parse_dates=True)

df['lifetime_err'] = get_errs(df.lifetime, df.anode_height, df.cathode_height, df.drift_time)
df['lifetime_smooth'] = medfilt(df.lifetime, kernel_size=501)
# clean
df_clean = df[
    (df.index < dt.datetime(2021, 1, 31, tzinfo=dt.timezone.utc))
    & (df.lifetime > 0)
]

df_clean['tdiff_days'] = (df_clean.index - df_clean.index[0]).total_seconds()/60./60./24.


def get_RHS(x, t, p):
    rhs = p[1]/(1+(t/20))
    rhs += -1.0*x/p[2]
    rhs += -1.0*x/p[0]
    return rhs


outgassing = inv_outg(0.075e-3)  # g/d
flow_speed = 60.0
lxe_flow_speed = 2.0  # lpm
efficiency = 0.57

tmax = 12.0
x0 = 0.00062


times = np.linspace(0, tmax, 10000)

tau = lxe_mass/density*1000/flow_speed*60/3600/24  # d
tau2 = lxe_mass/lxe_density/lxe_flow_speed/efficiency*60/3600/24  # d
print('outgassing', 'tau', 'tau2')
print(outgassing, tau, tau2)
p = [tau, outgassing / (lxe_mass/lxe_density), tau2]
sol = odeint(
    get_RHS,
    x0,
    times,
    args=(p,)
)
p = [tau, outgassing / (lxe_mass/lxe_density), tau2]
sol = odeint(
    get_RHS,
    x0,
    times,
    args=(p,)
)
fig = plt.figure(figsize=(12, 8))
plt.plot(df_clean.tdiff_days, df_clean.lifetime_smooth, 'k.', rasterized=True)
plt.fill_between(df_clean.tdiff_days, df_clean.lifetime_smooth-df_clean.lifetime_err, df_clean.lifetime_smooth+df_clean.lifetime_err, color='slategrey', alpha=0.3)
plt.plot(times, 1.0/sol.T[0], 'C0--', linewidth=3)
plt.xlim([0, 11])
plt.ylim([1000, 8000])
#plt.yscale('log')
plt.grid(alpha=0.3)
plt.xlabel('Time [d]')
plt.ylabel('Electron Lifetime [$\mu$s]')
if True:
    plt.savefig('plots/st707_efficiency_trend_220526.pdf')
    plt.savefig('plots/st707_efficiency_trend_220526.png', dpi=200)
plt.show()

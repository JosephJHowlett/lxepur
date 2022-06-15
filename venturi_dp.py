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

from scipy.optimize import curve_fit

def string_to_lngs_epoch(string, datetime_format='%Y-%m-%d %H:%M:%S'):
    # Just putting things in LNGS time and making a unified date format
    # to make things simple.
    epoch = datetime.datetime(1970, 1, 1)
    return (datetime.datetime.strptime(string, datetime_format) - epoch).total_seconds()

flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)
pressure_offsets_old = {
    'pt101':0.0000, 'pt271':+0.0380, 'pt273':+0.0296, 'pt275':+0.0036, 'pt286':+0.0237, 'pt287':+0.0511
}
pressure_offsets = {
    'pt101':0.0000, 'pt271':+0.0400, 'pt275':+0.0052, 'pt286':+0.0153, 'pt287':+0.0559, 'pt289':+0.0689
}

def fix_pressures(df):
    for key, val in pressure_offsets.items():
        if key in df.columns.values:
            print('adjusting ', key)
            df[key] -= val
    return

def h_to_dP(h, rho=2818.0):
    g = 9.8
    return rho*g*h/1e5  # bar

def dP_to_h(dP, rho=2818.0):
    g = 9.8
    return (dP*1e5)/(rho*g)  # m


use_scada = False
if use_scada:
    from straxen import SCADAInterface
    sc = SCADAInterface()

    start_time = '2020-01-30 10:30:00'
    stop_time = '2020-02-09 10:30:00'

    var_names = {
        'venturi_dp': 'XE1T.PUR_DPT288_PMON_A.PI',
        'venturi_flow': 'XE1T.PUR_DPT288_FLOW_A.PI',
        'pt286': 'XE1T.PUR_PT286_POUTLET.PI',
        'pt289': 'XE1T.PUR_PT289_PRADONCOLUMNINLET.PI',
        'te281': 'XE1T.PUR_TE281_TFEOUTLET.PI'
    }

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    fix_pressures(df)
    df['venturi_deltaP'] = df['pt286'] - df['pt289']
    df['venturi_flow'] = flow_func(df['venturi_dp'])
    df.to_csv('data/2020-01-30_2020-02-09_venturi_data-fixed.csv')

    ## long period with unstable, high flow

    start_time = '2021-02-01 10:30:00'
    stop_time = '2021-03-18 10:30:00'

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    fix_pressures(df)
    df['venturi_deltaP'] = df['pt286'] - df['pt289']
    df['venturi_flow'] = flow_func(df['venturi_dp'])
    df.to_csv('data/2021-02-01_2021-03-18_venturi_data-fixed.csv')

    ## long period stable 2.0 lpm

    start_time = '2021-04-15 10:30:00'
    stop_time = '2022-01-02 10:30:00'

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=1000,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    fix_pressures(df)
    df['venturi_deltaP'] = df['pt286'] - df['pt289']
    df['venturi_flow'] = flow_func(df['venturi_dp'])
    df.to_csv('data/2021-04-15_2022-01-02_venturi_data-fixed.csv')

    ## commissioning

    start_time = '2020-10-13 10:30:00'
    stop_time = '2020-12-31 10:30:00'

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=1000,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    fix_pressures(df)
    df['venturi_deltaP'] = df['pt286'] - df['pt289']
    df['venturi_flow'] = flow_func(df['venturi_dp'])
    df.to_csv('data/2020-10-13_2020-12-31_venturi_data-fixed.csv')


rho = 2.88
h = 1.7
pmin = 0.3
pmin = 0
pmin = -0.05

offset = -1*-0.03835128423822784


df = pd.read_csv('data/2021-02-01_2021-03-18_venturi_data-fixed.csv')
df['venturi_deltaP'] -= h_to_dP(h)
df_cut = df[
    #(df['venturi_deltaP'] > 0.3)
    (df['venturi_flow'] > 0.5)
    & (df['venturi_deltaP'] > pmin)
]
#plt.scatter(df_cut['venturi_flow'], df_cut['venturi_deltaP'] + offset, s=3, color='C1', label='Unstable flow period (March 2021)')

df = pd.read_csv('data/2021-04-15_2022-01-02_venturi_data-fixed.csv')
df['venturi_deltaP'] -= h_to_dP(h)
df_cut = df[
    #(df['venturi_deltaP'] > 0.3)
    (df['venturi_flow'] > 0.5)
    & (df['venturi_deltaP'] > pmin)
]
#plt.scatter(df_cut['venturi_flow'], df_cut['venturi_deltaP'] + offset, s=3, color='C2', label='Steady operations (2021-2022)')

print(df_cut['venturi_flow'].mean(), df_cut['venturi_deltaP'].mean())

df = pd.read_csv('data/2020-01-30_2020-02-09_venturi_data-fixed.csv')
df['venturi_deltaP'] -= h_to_dP(h)
df_cut = df[
    (df['venturi_flow'] > 0.5)
    & (df['venturi_deltaP'] > pmin)
]



flowcoef = lambda x, a, c: a*x**2.0 + c
#flowcoef = lambda x, a: a*x**2.0

popt, popcov = curve_fit(flowcoef, df_cut['venturi_flow'], df_cut['venturi_deltaP'])

plt.scatter(df_cut['venturi_flow'], df_cut['venturi_deltaP'] + offset, s=3, color='C0', label='Second commissioning test (Feb 2020)')
xx = np.linspace(df_cut['venturi_flow'].values.min(), df_cut['venturi_flow'].values.max(), 1000)
plt.plot(xx, flowcoef(xx, *popt) + offset, 'C0--')
print(popt)

df = pd.read_csv('data/2020-10-13_2020-12-31_venturi_data-fixed.csv')
df['venturi_deltaP'] -= h_to_dP(h)
df_cut = df[
    #(df['venturi_deltaP'] > 0.3)
    (df['venturi_flow'] > 0.5)
    & (df['venturi_deltaP'] > pmin)
]
#plt.scatter(df_cut['venturi_flow'], df_cut['venturi_deltaP'] + offset, s=3, color='C3', label='Commissioning (Oct-Dec 2020)')
plt.ylabel('Venturi Pressure Drop [bar]')
plt.xlabel('Venturi Flow [liter / min]')
#plt.legend(loc=2)
plt.ylim([0, 0.15])
if True:
    plt.savefig('plots/venturi_dp.png', dpi=300)
    plt.savefig('plots/venturi_dp.pdf')
plt.show()

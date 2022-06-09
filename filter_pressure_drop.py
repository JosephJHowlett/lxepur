import matplotlib.pyplot as plt
import matplotlib


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

import matplotlib.dates as mdates
myFmt = mdates.DateFormatter('%Y-%m-%d %H:%M:%S')

flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)
pressure_offsets = {
    'pt101':0.0000, 
    'pt271':+0.0380, 
    'pt273':+0.0296, 
    'pt275':+0.0036, 
    'pt286':+0.0237, 
    'pt287':+0.0511
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
data_filename = 'data/data_filter_dp_CP2.csv'

if use_scada:
    from straxen import SCADAInterface
    sc = SCADAInterface()

    start_time = '2020-01-30 10:30:00'
    stop_time = '2020-02-09 10:30:00'

    var_names = {
        'venturi_dp': 'XE1T.PUR_DPT288_PMON_A.PI',
        'venturi_flow': 'XE1T.PUR_DPT288_FLOW_A.PI',
        'pt286': 'XE1T.PUR_PT286_POUTLET.PI',
        'pt274': 'XE1T.PUR_PT274_POXYGENFILTERINLET.PI',
        'pt273': 'XE1T.PUR_PT273_POXYTGENFILTERINLET.PI',
        'fv284': 'XE1T.PUR_FV284_B_DOUT_.ST',
        'fv283': 'XE1T.PUR_FV283_B_DOUT_.ST',
        'fv275': 'XE1T.PUR_FV275_B_DOUT_.ST',
        'fv282': 'XE1T.PUR_FV282_B_DOUT_.ST',
        'fv277': 'XE1T.PUR_FV277_B_DOUT_.ST',
        'fv278': 'XE1T.PUR_FV278_B_DOUT_.ST',
        'fv279': 'XE1T.PUR_FV279_B_DOUT_.ST',
        'fv280': 'XE1T.PUR_FV280_B_DOUT_.ST',
        #'pt289': 'XE1T.PUR_PT289_PRADONCOLUMNINLET.PI',
        
    }

    df = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    fix_pressures(df)
    df['venturi_flow'] = flow_func(df['venturi_dp'])

    ## tube v2
    start_time = '2020-11-06 07:43:00'
    stop_time = '2020-11-06 15:30:00'

    var_names = {
        'venturi_dp': 'XE1T.PUR_DPT288_PMON_A.PI',
        'venturi_flow': 'XE1T.PUR_DPT288_FLOW_A.PI',
        'pt286': 'XE1T.PUR_PT286_POUTLET.PI',
        'pt274': 'XE1T.PUR_PT274_POXYGENFILTERINLET.PI',
        'pt273': 'XE1T.PUR_PT273_POXYTGENFILTERINLET.PI',
        'fv284': 'XE1T.PUR_FV284_B_DOUT_.ST',
        'fv283': 'XE1T.PUR_FV283_B_DOUT_.ST',
        'fv275': 'XE1T.PUR_FV275_B_DOUT_.ST',
        'fv282': 'XE1T.PUR_FV282_B_DOUT_.ST',
        'fv277': 'XE1T.PUR_FV277_B_DOUT_.ST',
        'fv278': 'XE1T.PUR_FV278_B_DOUT_.ST',
        'fv279': 'XE1T.PUR_FV279_B_DOUT_.ST',
        'fv280': 'XE1T.PUR_FV280_B_DOUT_.ST',
        #'pt289': 'XE1T.PUR_PT289_PRADONCOLUMNINLET.PI',
        
    }

    df_tube_v2 = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=100,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    
    fix_pressures(df_tube_v2)
    df_tube_v2['venturi_flow'] = flow_func(df_tube_v2['venturi_dp'])

    ## St 707 v1 (SR0)
    start_time = '2020-12-09 00:00:00'
    stop_time = '2020-12-10 00:00:00'

    var_names = {
        'venturi_dp': 'XE1T.PUR_DPT288_PMON_A.PI',
        'venturi_flow': 'XE1T.PUR_DPT288_FLOW_A.PI',
        'pt286': 'XE1T.PUR_PT286_POUTLET.PI',
        'pt274': 'XE1T.PUR_PT274_POXYGENFILTERINLET.PI',
        'pt273': 'XE1T.PUR_PT273_POXYTGENFILTERINLET.PI',
        'fv284': 'XE1T.PUR_FV284_B_DOUT_.ST',
        'fv283': 'XE1T.PUR_FV283_B_DOUT_.ST',
        'fv275': 'XE1T.PUR_FV275_B_DOUT_.ST',
        'fv282': 'XE1T.PUR_FV282_B_DOUT_.ST',
        'fv277': 'XE1T.PUR_FV277_B_DOUT_.ST',
        'fv278': 'XE1T.PUR_FV278_B_DOUT_.ST',
        'fv279': 'XE1T.PUR_FV279_B_DOUT_.ST',
        'fv280': 'XE1T.PUR_FV280_B_DOUT_.ST',
        #'pt289': 'XE1T.PUR_PT289_PRADONCOLUMNINLET.PI',
        
    }

    df_st707_sr0 = sc.get_scada_values(
        var_names,
        start=string_to_lngs_epoch(start_time)*1e9,
        end=string_to_lngs_epoch(stop_time)*1e9,
        every_nth_value=10,
        fill_gaps="forwardfill", # the DB stops recording values that do not change
        query_type_lab=True,
    )
    
    fix_pressures(df_st707_sr0)
    df_st707_sr0['venturi_flow'] = flow_func(df_st707_sr0['venturi_dp'])
else:
    # get data from csvs
    df = pd.read_csv('data/data_filter_dp_CP2.csv', index_col='time UTC', parse_dates=True)
    df_tube_v2 = pd.read_csv('data/data_filter_dp_tube_v2.csv', index_col='time UTC', parse_dates=True)
    df_st707_sr0 = pd.read_csv('data/data_filter_dp_st707_v1.csv', index_col='time UTC', parse_dates=True)



valve_states = {}
for col in df.columns:
    if 'fv' in col:
        valve_states[col] = df[col].values.astype(int)

tpc_mode = (
    valve_states['fv275'] 
    & valve_states['fv284']
    & (~ valve_states['fv282'])
    & (~ valve_states['fv283'])
)
filter_mode = (
    (~ valve_states['fv275'] )
    & (~valve_states['fv284'])
    & valve_states['fv282']
    & valve_states['fv283']
)
filter_271 = (
    valve_states['fv279'] 
    & valve_states['fv280']
)
filter_274 = (
    valve_states['fv277'] 
    & valve_states['fv278']
)

df_tpc_mode_271 = df[
    (tpc_mode==1) & (df['venturi_flow']>0)
    & (filter_271==1)
]
df_tpc_mode_274 = df[
    (tpc_mode==1) & (df['venturi_flow']>0)
    & (filter_274==1)
]
df_filter_mode_271 = df[
    (filter_mode==1) 
    & (df['venturi_flow']>0)
    & (filter_271==1)
]
df_filter_mode_274 = df[
    (filter_mode==1) 
    & (df['venturi_flow']>0)
    & (filter_274==1)
]

quadfit = lambda x, a, b: a*(x**2.0) + b*x
quadfit_c = lambda x, a, b, c: a*(x**2.0) + b*x + c
h = 0.258  # m
early_time = lambda df: df[df.index < '2020-02-07 10:43:20+00:00']
late_time = lambda df: df[df.index > '2020-02-07 10:43:20+00:00']

## pancake v2
df_pancake_v2 = pd.concat([early_time(df_filter_mode_274), early_time(df_tpc_mode_274)])
df_pancake_v2 = df_pancake_v2[(df_pancake_v2['venturi_flow'] > 1.0) | (df_pancake_v2['pt273'] - df_pancake_v2['pt286'] + h_to_dP(h) < 0.2)]

## st707 v0
df_st707_v0 = pd.concat([early_time(df_filter_mode_271), early_time(df_tpc_mode_271)])
df_st707_v0 = df_st707_v0[df_st707_v0['venturi_flow'] < 3]
# remove outliers
popt, popcov = curve_fit(quadfit_c, df_st707_v0['venturi_flow'], df_st707_v0['pt274'] - df_st707_v0['pt286'] + h_to_dP(h))
df_st707_v0 = df_st707_v0[
    np.absolute(df_st707_v0['pt274'] - df_st707_v0['pt286'] + h_to_dP(h) - quadfit_c(df_st707_v0['venturi_flow'], *popt))/quadfit_c(df_st707_v0['venturi_flow'], *popt)<0.25
]

## tube v1
df_tube_v1 = pd.concat([late_time(df_filter_mode_274), late_time(df_tpc_mode_274)])
popt, popcov = curve_fit(quadfit, df_tube_v1['venturi_flow'], df_tube_v1['pt273'] - df_tube_v1['pt286'] + h_to_dP(h))
df_tube_v1 = df_tube_v1[np.absolute(df_tube_v1['pt273'] - df_tube_v1['pt286'] + h_to_dP(h) - quadfit(df_tube_v1['venturi_flow'], *popt))<0.1]


## cleanup st707 v1
popt, popcov = curve_fit(quadfit, df_st707_sr0['venturi_flow'], df_st707_sr0['pt274'] - df_st707_sr0['pt286'] + h_to_dP(h))

df_st707_sr0 = df_st707_sr0[quadfit(df_st707_sr0['venturi_flow'], *popt) - (df_st707_sr0['pt274'] - df_st707_sr0['pt286'] + h_to_dP(h)) < 0.07]
popt, popcov = curve_fit(
    quadfit, df_st707_sr0['venturi_flow'], df_st707_sr0['pt274'] - df_st707_sr0['pt286'] + h_to_dP(h)
)


## redefine fits for useful parameter units
quadfit = lambda x, a, b: a*(x**2.0) + b*x
quadfit = lambda x, a, b: (x**2.0)/(a**2.0) + x/b

# pancake v2
plt.scatter(
    df_pancake_v2['venturi_flow'], 
    df_pancake_v2['pt273'] - df_pancake_v2['pt286'] + h_to_dP(h), 
)
popt, popcov = curve_fit(quadfit, df_pancake_v2['venturi_flow'], df_pancake_v2['pt273'] - df_pancake_v2['pt286'] + h_to_dP(h))
xs = np.linspace(0, 3.75, 1000)
plt.plot(xs, quadfit(xs, *popt), label='Filter C1')
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt, np.sqrt(np.diag(popcov)))])

# st707 v0
plt.scatter(
    df_st707_v0['venturi_flow'], 
    df_st707_v0['pt274'] - df_st707_v0['pt286'] + h_to_dP(h), 
)
popt, popcov = curve_fit(quadfit, df_st707_v0['venturi_flow'], df_st707_v0['pt274'] - df_st707_v0['pt286'] + h_to_dP(h))
xs = np.linspace(0, 2.5, 1000)
plt.plot(xs, quadfit(xs, *popt), label='Filter C2')
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt, np.sqrt(np.diag(popcov)))])

# tube v1
popt, popcov = curve_fit(quadfit, df_tube_v1['venturi_flow'], df_tube_v1['pt273'] - df_tube_v1['pt286'] + h_to_dP(h))
plt.scatter(
    df_tube_v1['venturi_flow'], 
    df_tube_v1['pt273'] - df_tube_v1['pt286'] + h_to_dP(h), 
)

xs = np.linspace(0, 2.5, 1000)
plt.plot(xs, quadfit(xs, *popt), label='Filter C3')
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt, np.sqrt(np.diag(popcov)))])

# tube v2

plt.scatter(
    df_tube_v2['venturi_flow'], 
    df_tube_v2['pt273'] - df_tube_v2['pt286'] + h_to_dP(h), 
)
popt, popcov = curve_fit(quadfit, df_tube_v2['venturi_flow'], df_tube_v2['pt273'] - df_tube_v2['pt286'] + h_to_dP(h))
xs = np.linspace(-0.5, 3.75, 1000)
plt.plot(xs, quadfit(xs, *popt), label='Filter 1')
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt, np.sqrt(np.diag(popcov)))])

# st707 sr0
popt, popcov = curve_fit(
    quadfit, df_st707_sr0['venturi_flow'], df_st707_sr0['pt274'] - df_st707_sr0['pt286'] + h_to_dP(h)
)
plt.scatter(
    df_st707_sr0['venturi_flow'], 
    df_st707_sr0['pt274'] - df_st707_sr0['pt286'] + h_to_dP(h), 
)
xs = np.linspace(-0.5, 3.75, 1000)
plt.plot(xs, quadfit(xs, *popt), label='Filter 2')
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt, np.sqrt(np.diag(popcov)))])

plt.xlabel('LXe Flow [liter / min]')
plt.ylabel('Pressure Drop [bar]')
plt.grid(alpha=0.3)

plt.xlim([0, 4])
plt.ylim([0, 1])
plt.legend()
if True:
    plt.savefig('plots/filter_dp.png', dpi=300)
    plt.savefig('plots/filter_dp.pdf')
plt.show()

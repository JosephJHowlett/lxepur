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

import time, calendar
from scipy.optimize import curve_fit

def string_to_lngs_epoch(string, datetime_format='%m-%d-%y %H:%M:%S'):
    # Just putting things in LNGS time and making a unified date format
    # to make things simple.
    epoch = datetime.datetime(1970, 1, 1)
    return (datetime.datetime.strptime(string, datetime_format) - epoch).total_seconds()

flow_func = lambda dp: (-0.8501+np.sqrt(0.8501**2+4.*22.573*dp))/(2.*22.573)

def convert_fitval(fitval):
    Cp = 346.24  # J/kg/C
    rho = 2.8188  # kg/l
    return fitval*rho*Cp/60  # W

def convert_inverse_fitval(fitval):
    Cp = 346.24  # J/kg/C
    rho = 2.8188  # kg/l
    return (rho*Cp/60)/fitval  # W

def get_means_stds_from_viewercsv(csv_files, time_intervals):
    to_unixtime = lambda t: calendar.timegm(time.strptime(t[:17], '%m/%d/%y %H:%M:%S'))-1*3600
    dataframes_list = [pd.read_csv(f, converters={'LNGS Time': lambda t: calendar.timegm(time.strptime(t, '%Y-%m-%d %H:%M:%S'))-1*3600}) for f in csv_files]
    sensors_dataframe = pd.concat(dataframes_list)
    unixtime_vals_array = sensors_dataframe['LNGS Time'].values.astype(np.float64)
    df = []
    for time_interval in time_intervals:
        t0, t1 = tuple(map(to_unixtime, time_interval))
        interval_array = (unixtime_vals_array >= t0) & (unixtime_vals_array <= t1)
        row = {}
        for col in sensors_dataframe.columns:
            if 'LNGS' in col:
                continue
            colarr = sensors_dataframe[col].values.astype(np.float64)[interval_array]
            row[col + '_mean'] = colarr.mean()
            row[col + '_std'] = colarr.std()
        df.append(row)
    df = pd.DataFrame(df)
    return df


fitfunc = lambda flow, a, c: a/flow + c

## filter mode data
dpt288 = np.array([52.03, 73.65, 50.13, 99.92, 23.99, 50.89])
flow_274_281 = flow_func(dpt288)
te274 = np.array([-90.031, -90.001, -89.494, -90.550, -88.884, -89.928])
te281 = np.array([-89.885, -89.897, -89.323, -90.556, -88.322, -89.751])
tdiff_274_281 = te281 - te274
popt_274_281, popcov = curve_fit(
    fitfunc, 
    flow_274_281, 
    tdiff_274_281, 
)
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt_274_281, np.sqrt(np.diag(popcov)))])
xs = np.linspace(0.8, 2.2, 1000)
plt.plot(xs, fitfunc(xs, *popt_274_281), 'k--')
plt.plot(flow_274_281, tdiff_274_281, 'k.')
plt.xlabel('venturi flow [liter / min]')
plt.ylabel(r'$T_{281}^m - T_{274}^m}$')
print(convert_fitval(popt_274_281[0]))
plt.show()

## filter 271 in TPC mode
csv_files = [
    'data/xenonnt-12_06_2021-12_28_2021.csv',
    'data/xenonnt-03_09_2021-03_29_2021.csv',
    'data/xenonnt-04_21_2021-04_30_2021.csv',
    'data/xenonnt-12_01_2020-01_30_2021.csv',
]
time_intervals = [
    ['12/10/21 10:11:00', '12/11/21 06:23:00'],
    ['12/14/21 10:46:00', '12/15/21 05:30:00'],
    #['12/16/21 12:03:00', '12/16/21 16:29:00'],
    #['12/19/21 21:21:00', '12/20/21 07:42:00'],
    ['12/21/21 06:51:00', '12/21/21 12:46:00'],
    ['12/24/21 05:49:00', '12/24/21 14:41:00'],
    ['12/27/21 15:37:00', '12/28/21 01:28:00'],
 
    ['03/14/21 10:05:00', '03/19/21 02:05:00'],
    ['03/22/21 13:05:00', '03/23/21 11:05:00'],
    ['03/28/21 07:05:00', '03/29/21 05:05:00'],
    
    ['04/28/21 07:05:00', '04/29/21 05:05:00'],
    
    ['12/07/20 03:05:00', '12/09/20 04:05:00'],
    ['12/26/20 03:05:00', '12/27/20 04:05:00'],
    ['01/05/21 11:05:00', '01/07/21 03:05:00'],
    #['01/20/21 03:05:00', '01/20/21 09:05:00'],
]
combined_df = get_means_stds_from_viewercsv(csv_files, time_intervals)
flow_278_281 = flow_func(combined_df['PUR_DPT288_PMON_A.PI_mean'])
tdiff_278_281 = combined_df['PUR_TE281_TFEOUTLET.PI_mean'] - combined_df['PUR_TE278_TFILTEROUTLET.PI_mean']

popt_278_281, popcov = curve_fit(
    fitfunc, 
    flow_278_281, 
    tdiff_278_281, 
)
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt_278_281, np.sqrt(np.diag(popcov)))])
xs = np.linspace(1.2, 3, 1000)
plt.plot(xs, fitfunc(xs, *popt_278_281), 'k--')

plt.plot(
    flow_278_281,
    tdiff_278_281,
    'k.'
)
plt.xlabel('venturi flow [liter / min]')
plt.ylabel(r'$T_{281}^m - T_{278}^m}$')
print(convert_fitval(popt_278_281[0]))
plt.show()

## filter 274 in TPC mode
csv_files = [
    'data/xenonnt-11_07_2020-11_30_2020.csv',
    'data/xenonnt-04_17_2021-04_21_2021.csv',
    'data/xenonnt-02_22_2022-03_21_2022.csv',
]
time_intervals = [
    ['11/14/20 18:05:00', '11/16/20 04:05:00'],
    ['11/28/20 04:05:00', '11/29/20 18:05:00'],
    ['04/20/21 00:05:00', '04/20/21 22:05:00'],
    
    ['02/23/22 00:06:00', '02/23/22 19:05:00'],
    #['03/02/22 18:06:00', '03/03/22 11:05:00'],
    #['03/15/22 09:06:00', '03/16/22 10:05:00'],
    #['03/19/22 01:06:00', '03/19/22 16:05:00'],
    #['03/02/22 18:06:00', '03/20/22 11:05:00'],
]
combined_df = get_means_stds_from_viewercsv(csv_files, time_intervals)

flow_276_281 = flow_func(combined_df['PUR_DPT288_PMON_A.PI_mean'])
tdiff_276_281 = combined_df['PUR_TE281_TFEOUTLET.PI_mean'] - combined_df['PUR_TE276_TFILTEROUTLET.PI_mean']

popt_276_281, popcov = curve_fit(
    fitfunc, 
    flow_276_281, 
    tdiff_276_281, 
)
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt_276_281, np.sqrt(np.diag(popcov)))])
xs = np.linspace(1.5, 3.5, 1000)
plt.plot(xs, fitfunc(xs, *popt_276_281), 'k--')

plt.plot(
    flow_276_281,
    tdiff_276_281,
    'k.'
    #xerr=flow_func(combined_df['PUR_DPT288_PMON_A.PI_std']),
    #yerr=np.sqrt(combined_df['PUR_TE281_TFEOUTLET.PI_std']**2.0 + combined_df['PUR_TE278_TFILTEROUTLET.PI_std']**2.0),
    #fmt='k.'
)
plt.xlabel('venturi flow [liter / min]')
plt.ylabel(r'$T_{281}^m - T_{276}^m}$')
print(convert_fitval(popt_276_281[0]))
plt.show()


## subtract offsets, combine
L_274_281 = 2.743  # m
L_276_281 = 3.6  # m
L_278_281 = 3.6  # m

plt.plot(
    flow_274_281,
    L_274_281/(tdiff_274_281 - popt_274_281[1]),
    '.',
    label=r'TE-274 $\rightarrow$ TE-281',
    #xerr=flow_func(combined_df['PUR_DPT288_PMON_A.PI_std']),
    #yerr=np.sqrt(combined_df['PUR_TE281_TFEOUTLET.PI_std']**2.0 + combined_df['PUR_TE278_TFILTEROUTLET.PI_std']**2.0),
    #fmt='k.'
)
plt.plot(
    flow_278_281,
    L_278_281/(tdiff_278_281 - popt_278_281[1]),
    '.',
    label=r'TE-278 $\rightarrow$ TE-281',
    #xerr=flow_func(combined_df['PUR_DPT288_PMON_A.PI_std']),
    #yerr=np.sqrt(combined_df['PUR_TE281_TFEOUTLET.PI_std']**2.0 + combined_df['PUR_TE278_TFILTEROUTLET.PI_std']**2.0),
    #fmt='k.'
)
plt.plot(
    flow_276_281,
    L_276_281/(tdiff_276_281 - popt_276_281[1]),
    '.',
    label=r'TE-276 $\rightarrow$ TE-281',
    #xerr=flow_func(combined_df['PUR_DPT288_PMON_A.PI_std']),
    #yerr=np.sqrt(combined_df['PUR_TE281_TFEOUTLET.PI_std']**2.0 + combined_df['PUR_TE278_TFILTEROUTLET.PI_std']**2.0),
    #fmt='k.'
)

line = lambda flow, a: a*flow

popt, popcov = curve_fit(
    line,
    flow_274_281.tolist() + flow_278_281.tolist() + flow_276_281.tolist(),
    (L_274_281/(tdiff_274_281 - popt_274_281[1])).tolist() + (L_278_281/(tdiff_278_281 - popt_278_281[1])).tolist() + (L_276_281/(tdiff_276_281 - popt_276_281[1])).tolist(),
)
xs = np.linspace(1, 3, 100)
plt.plot(xs, line(xs, *popt), 'k--')
print(['%.4f ± %.4f' % (val, err) for val, err in zip(popt, np.sqrt(np.diag(popcov)))])
plt.legend()
plt.ylabel(r'$\Delta L/\Delta T_\mathrm{true}$ [m/K]')
plt.xlabel('Venturi Flow [l/min]')
plt.grid(alpha=0.3)
if True:
    plt.savefig('plots/dqdl_all_combined.pdf')
plt.show()

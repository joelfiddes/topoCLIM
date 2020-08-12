import matplotlib
import matplotlib.pyplot as plt
import melodist
import numpy as np
import pandas as pd
import scipy.stats

#need tmin and max!!!

path_inp = "/home/joel/sim/qmap/wfj_long.csv"
df= pd.read_csv(path_inp, index_col=0, parse_dates=True)

# resample 3h OBS to 1h data (psum no longer valid) - can use 1h data directly here in future
df_interpol = df.resample('H').mean()
df_interpol = df_interpol.interpolate()
df_interpol.head(4)

    
df2 = df_interpol[['TA','PINT', 'RH', 'ISWR', 'VW']]  
df2.rename(columns={'TA': 'temp', 'PINT': 'precip', 'RH': 'hum', 'ISWR': 'glob', 'VW': 'wind'}, inplace=True)

# unit conversions
df2.precip *=1#kgm-2s-1 to mm/hr
df2.hum *=100

def plot(obs, sim):
    plt.figure()
    ax = plt.gca()
    obs.loc[plot_period].plot(ax=ax, color='black', label='obs', lw=2)
    sim.loc[plot_period].plot(ax=ax)
    plt.legend()
    plt.show()

def calc_stats(obs, sim):
    df = pd.DataFrame(columns=['mean', 'std', 'r', 'rmse', 'nse'])

    obs = obs.loc[validation_period]
    sim = sim.loc[validation_period]
    
    df.loc['obs'] = obs.mean(), obs.std(), 1, 0, 1
    
    for c in sim.columns:
        osdf = pd.DataFrame(data=dict(obs=obs, sim=sim[c])).dropna(how='any')
        o = osdf.obs
        s = osdf.sim
    
        r = scipy.stats.pearsonr(o, s)[0]
        rmse = np.mean((o - s)**2)
        nse = 1 - np.sum((o - s)**2) / np.sum((o - o.mean())**2)
        df.loc[c] = s.mean(), s.std(), r, rmse, nse
    
    return df



def print_stats(obs, sim):
    df = calc_stats(obs, sim)
    html = df.round(2).style
    return html


# station attributes
longitude = 9.80988
latitude = 46.82786
timezone = 1



calibration_period = slice('2000-01-01', '2015-12-31')
validation_period = slice('2016-01-01', '2016-12-31')
plot_period = slice('2016-09-03', '2016-09-05')
data_obs_hourly = df2
data_obs_hourly = melodist.util.drop_incomplete_days(data_obs_hourly)
data_obs_daily = melodist.util.daily_from_hourly(data_obs_hourly)


import xarray as xr
nc ='/home/joel/sim/qmap/test/pyout/aresult/NOAA-GFDL-GFDL-ESM2M_rcp85_r1i1p1_SMHI-RCA4_v1__TS.nc_TS_ALL_ll.nc'
ds = xr.open_dataset(nc,decode_times=True)
datetimeindex = ds.indexes["time"].to_datetimeindex() 
datetimeindex.is_all_dates  
#ds['dtime'] = datetimeindex
# myx =np.argmin(abs(ds['lon']-longitude))   
# myy = np.argmin(abs(ds['lat']-latitude))   

df = ds.sel(lon=longitude, lat=latitude, method='nearest').to_dataframe()



#df = ds.to_dataframe()
df2 = df[['tas','pr', 'hurs', 'rsds', 'uas']]#, 'vas']]  
df2.rename(columns={'tas': 'temp', 'pr': 'precip', 'hurs': 'hum', 'rsds': 'glob', 'uas': 'wind'}, inplace=True)


#df2["datetime"] = datetimeindex

# set index
df3 =df2.set_index(datetimeindex) 

# resample to daily (adds missing days i think - converts implicitly to standard calende not noleap)
df4 = df3.resample('D').mean()


station = melodist.Station(lon=longitude, lat=latitude, timezone=timezone, data_daily=df4)

# aggregated hourly obs with min and max for calculating stats
station2 = melodist.Station(lon=longitude, lat=latitude, timezone=timezone, data_daily=data_obs_daily)


station2.statistics = melodist.StationStatistics(data_obs_hourly.loc[calibration_period])

stats = station2.statistics
stats.calc_wind_stats()
stats.calc_humidity_stats()
stats.calc_temperature_stats()
stats.calc_radiation_stats()
stats.calc_precipitation_stats()


tempdf = pd.DataFrame()
for method in ('sine_min_max', 'sine_mean', 'mean_course_min_max', 'mean_course_mean'):
    station.disaggregate_temperature(method=method, min_max_time='sun_loc_shift')
    tempdf[method] = station.data_disagg.temp


plot(data_obs_hourly.temp, tempdf)
print_stats(data_obs_hourly.temp, tempdf)

# ok, let's say we have decided to use sine_min_max. next we want to disaggregate humidity.
# as some of hum disagg functions rely on disagg'd temperature values we disagg temp again
# with our chosen method
station2.disaggregate_temperature(method='mean_course_mean', min_max_time='sun_loc_shift')


humdf = pd.DataFrame()
for method in ('equal', 'minimal', 'dewpoint_regression',
               'linear_dewpoint_variation', 'min_max',
               'month_hour_precip_mean'):
    station.disaggregate_humidity(method=method)
    humdf[method] = station.data_disagg.hum

plot(data_obs_hourly.hum, humdf)
print_stats(data_obs_hourly.hum, humdf)

globdf = pd.DataFrame()
for method in ('pot_rad',
               # 'pot_rad_via_ssd',  # not possible here as we do not have sunshine duration data
               'pot_rad_via_bc',
               'mean_course'):
    station.disaggregate_radiation(method=method)
    globdf[method] = station.data_disagg.glob

plot(data_obs_hourly.glob, globdf)
print_stats(data_obs_hourly.glob, globdf)


winddf = pd.DataFrame()
for method in ('equal', 'cosine', 'random'):
    station.disaggregate_wind(method=method)
    winddf[method] = station.data_disagg.wind

plot(data_obs_hourly.wind, winddf)
print_stats(data_obs_hourly.wind, winddf)


precipdf = pd.DataFrame()
for method in ('equal', 'cascade'):
    station.disaggregate_precipitation(method=method)
    precipdf[method] = station.data_disagg.precip
    
plot(data_obs_hourly.precip, precipdf)




#final
station.disaggregate_temperature(method='sine_mean')
station.disaggregate_humidity(method='min_max')
station.disaggregate_radiation(method='mean_course')
station.disaggregate_wind(method='random')
station.disaggregate_precipitation(method='cascade')

# make wind direction equal
station.data_disagg['DW'] = df_interpol.DW.resample('H').mean().interpolate(method='pad') 

# disagg lwin with T
df_interpol.ILWR

plot(data_obs_hourly.temp,data_obs_hourly.glob)

df_interpol_day = df.resample('D').mean()

# compute a daily scaling factor and apply to each hour
scalingFact = df_interpol_day.ILWR/df_interpol_day.TA
scalingFact_day = scalingFact.resample('H').mean()
scalingFact_day.interpolate(method='pad')   
station.data_disagg['ILWR']  = station.data_disagg.temp * scalingFact_day
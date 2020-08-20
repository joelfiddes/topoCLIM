import matplotlib.pyplot as plt
import melodist
import numpy as np
import pandas as pd
import scipy.stats
import glob
import matplotlib


#======================================================================================
# Settings
#======================================================================================

# Do we want to run all methods and plot to look at performance?
test_methods=False 

# station attributes
longitude = 9.80988
latitude = 46.82786
timezone = 1
slope  = 0

calibration_period = slice('2000-01-01', '2015-12-31')
validation_period = slice('2016-01-01', '2016-12-31')
plot_period = slice('2016-09-03', '2030-10-13')

#======================================================================================
# Observations (toposacle output)
# Units:
# DW degrees
# ILWR: wm-2
# ISWR: wm-2
# PINT: mm/hr
# PSUM: mm/timestep
# Rf: mm/hr
# Sf: mm/hr
# TA: K
# VW: m s-1
# RH: 0-1
# P: Pa
#======================================================================================
path_inp = "/home/joel/sim/qmap/wfj_long.csv"

# lesson of psum mess is dont use psum!!

df_obs= pd.read_csv(path_inp, index_col=0, parse_dates=True)

# coorect for toposcale error where PINT = PSUM at 3hr timestep. now it is mm/h 19/8/20
#psum is actually wrong! 20/8/20 eg this is too little:

# correction!!!!!!!!
df_obs.PSUM *= 3.

# resample 3h OBS to 1h data (psum no longer valid) - can use 1h data directly here in future
df_interpol = df_obs.resample('H').mean()
df_interpol = df_interpol.interpolate()
df_interpol.head(4)
df_interpol.PSUM *= 1/3. # devide by 3 as is psum over 3h
    


data_obs_hourly = df_interpol[['TA','PINT', 'RH', 'ISWR','ILWR', 'VW']]  
data_obs_hourly.rename(columns={'TA': 'temp', 'PINT': 'precip', 'RH': 'hum', 'ISWR': 'glob', 'ILWR': 'ilwr', 'VW': 'wind'}, inplace=True)


# unit conversions
data_obs_hourly.precip *=1.#(1/3600.)#mm/hr to to kgm-2s-1 
data_obs_hourly.hum *=100. # 0-1 to 0-100

#data_obs_hourly = df2
data_obs_hourly = melodist.util.drop_incomplete_days(data_obs_hourly)

"""Aggregates data (hourly to daily values) according to the characteristics
of each variable (e.g., average for temperature, sum for precipitation)"""
data_obs_daily = melodist.util.daily_from_hourly(data_obs_hourly) # 

#======================================================================================
# Methods
#======================================================================================
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




#======================================================================================
# Climate data to dissagg
# 'tas' K
# 'tasmin' K
# 'tasmax' K
# 'pr' kg m-2 s-1
# 'hurs' % 0-100
# 'rsds' W m-2
# 'rlds' W m-2
# 'uas' m s-1
# 'vas' m s-1
# 'ps' pa

#======================================================================================

import xarray as xr
nc ='/home/joel/sim/qmap/test/pyout/aresult/NOAA-GFDL-GFDL-ESM2M_rcp85_r1i1p1_SMHI-RCA4_v1__TS.nc_TS_ALL_ll.nc'
nc='/home/joel/sim/qmap/test/pyout/aresult/CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN53_v1__TS.nc_TS_ALL_ll.nc'

ncfiles = glob.glob('/home/joel/sim/qmap/test/pyout/aresult/'+ "*_TS_ALL_ll.nc")

for nc in ncfiles:

	ds = xr.open_dataset(nc,decode_times=True)

	# this bit depend on the model and calender i think 

	#either this:
	datetimeindex = ds.indexes["time"]#.to_datetimeindex() 
	datetimeindex.is_all_dates  

	#or this
	#datetimeindex = ds.indexes["time"].to_datetimeindex() 
	#datetimeindex.is_all_dates 

				#ds['dtime'] = datetimeindex
				# myx =np.argmin(abs(ds['lon']-longitude))   
				# myy = np.argmin(abs(ds['lat']-latitude))   

	df = ds.sel(lon=longitude, lat=latitude, method='nearest').to_dataframe()

	if df.shape[1]<14:
		#print "not enouigh vars"
		next   
	try:
	#df = ds.to_dataframe()
		df2 = df[['tas','tasmin' , 'tasmax','pr', 'hurs', 'rsds', 'rlds', 'ps']]#, 'vas']]  
		#print 'ok'
		print nc
		print df.columns
		break
	except:
		#print "no vars"
		#print df.columns
		next


#conmute wind speed and wind direction
df2['ws'] = np.sqrt(df.uas**2+df.vas**2)
df2['wd'] = (180 / np.pi) * np.arctan(df.uas/df.vas) + np.where(df.vas>0,180,np.where(df.uas>0,360,0))

# unit conversions
df2.pr *=3600*24 # kgm-2s-1 (same as mm /s)-> mm/day 
df2.rename(columns={'tas': 'temp','tasmax': 'tmax','tasmin': 'tmin', 'pr': 'precip', 'hurs': 'hum', 'rsds': 'glob', 'ws': 'wind'}, inplace=True)


#df2["datetime"] = datetimeindex

# set index
df3 =df2.set_index(datetimeindex) 


# resample to daily (adds missing days i think - converts implicitly to standard calende not noleap)
df4 = df3.resample('D').mean()


station = melodist.Station(lon=longitude, lat=latitude, timezone=timezone, data_daily=df4)

# aggregated hourly obs with min and max for calculating stats
#station2 = melodist.Station(lon=longitude, lat=latitude, timezone=timezone, data_daily=data_obs_daily)


station.statistics = melodist.StationStatistics(data_obs_hourly.loc[calibration_period])

stats = station.statistics
stats.calc_wind_stats()
stats.calc_humidity_stats()
stats.calc_temperature_stats()
stats.calc_radiation_stats()
stats.calc_precipitation_stats()

if test_methods is True:
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
	               'linear_dewpoint_variation', 
	               #'min_max', not poss as no min max in obs
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
	for method in ('equal'):
	    station.disaggregate_precipitation(method=method)
	    precipdf[method] = station.data_disagg.precip
	    
	plot(data_obs_hourly.precip, precipdf)




#final dissagregation with selected methods
station.disaggregate_temperature(method='sine_mean')
station.disaggregate_humidity(method='equal')
station.disaggregate_radiation(method='pot_rad_via_bc')
station.disaggregate_wind(method='random')
station.disaggregate_precipitation(method='equal')

# make wind direction equal each hour
station.data_disagg['DW'] = df4.wd.resample('H').mean().interpolate(method='pad') 

# make air pressure equal each hour
station.data_disagg['P'] = df4.ps.resample('H').mean().interpolate(method='pad') 



# disagg lwin with T

# compute a daily scaling factor and apply to each hour
scalingFact = df4.rlds/df4.temp
scalingFact_day = scalingFact.resample('H').mean()
scale_fact_day_interp = scalingFact_day.interpolate(method='pad')   
station.data_disagg['ILWR']  = station.data_disagg.temp * scale_fact_day_interp

#plot(data_obs_hourly.ilwr,station.data_disagg.ILWR)




#======================================================================================
# # partition prate to rain snow (mm/hr)
#======================================================================================
		
lowthresh=272.15
highthresh = 274.15
d = {'prate': station.data_disagg.precip, 'ta': station.data_disagg.temp}
df = pd.DataFrame(data=d)
snow = df.prate.where(df.ta<lowthresh) 
rain=df.prate.where(df.ta>=highthresh) 

mix1S = df.prate.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
mix1T = df.ta.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
mixSno=(highthresh - mix1T) / (highthresh-lowthresh)
mixRain=1-mixSno
addSnow=mix1S*mixSno
addRain=mix1S*mixRain

# nas to 0
snow[np.isnan(snow)] = 0 	
rain[np.isnan(rain)] = 0 
addRain[np.isnan(addRain)] = 0 
addSnow[np.isnan(addSnow)] = 0 

#======================================================================================
# # linearly reduce snow to zero in steep slopes
#if steepSnowReduce=="TRUE": # make this an option if need that in future
#======================================================================================

snowSMIN=30.
snowSMAX=80.
slope=slope

k= (snowSMAX-slope)/(snowSMAX - snowSMIN)

if slope<snowSMIN:
	k=1
if slope>snowSMAX:
	k=0

snowTot=(snow+addSnow) * k
rainTot=rain + addRain




#write out dataframe geotop



df_gtop= pd.DataFrame({	
"Prec":station.data_disagg.precip,
				"Ws":station.data_disagg.wind,
				"Wd":station.data_disagg.DW, 
						"RH":station.data_disagg.hum,#*0.01, #meteoio 0-1
	"Tair":station.data_disagg.temp -273.15, 
						"SW":station.data_disagg.glob, 
				"LW":station.data_disagg.ILWR, 

				
				})
df_gtop.index = df_gtop.index.strftime("%d/%m/%Y %H:00")   

df_gtop.index.name="Date"

	# fill outstanding nan in SW routine with 0 (night)

	# fileout=outDir+"/meteo"+"c"+str(i+1)+".csv"
df_gtop.to_csv(path_or_buf='/home/joel/sim/qmap/geotop_sim/meteo0001.txt',na_rep=-999,float_format='%.3f')
	# logging.info(fileout + " complete")





fsm_out= pd.DataFrame({	
 				"year": df_out.index.year, 
 				"month": df_out.index.month, 
				"day": df_out.index.day, 
				"hour": df_out.index.hour,
				"ISWR":station.data_disagg.glob, 
				"ILWR":station.data_disagg.ILWR, 
				"Sf":snowTot/(60*60), # prate in mm/hr to kgm2/s
				"Rf":rainTot/(60*60), # prate in mm/hr to kgm2/s
				"TA":station.data_disagg.temp, 
				"RH":station.data_disagg.hum,#*0.01, #meteoio 0-1
				"VW":station.data_disagg.wind,
				"P":station.data_disagg.P,
				
				
				})

fsm_out.to_csv(path_or_buf='/home/joel/sim/qmap/disagg.txt' ,na_rep=-999,float_format='%.3f', header=False, sep='\t', index=False, 
	columns=['year','month','day', 'hour', 'ISWR', 'ILWR', 'Sf', 'Rf', 'TA', 'RH', 'VW', 'P'])


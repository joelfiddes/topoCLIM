
"""
python3

This module takes post processed daily CORDEX downloads from esgf_post.py and 
produces point timeseries of dissagreated hourly files with standard calenders

Example:
        $ python topoCLIM.py

Vars:
	all parameters set in script currently

Details:

	1. select nearest gid cel tp POI
	2. check for full suit of variables
	3. convert calenders
	4, dissagregate daily to hourly

	- converyts all types of no standard (no leap 360 day) cal to standard cal 
		adapted methods from xclim package
	- daily to hourly dissagregation, adapted methods from melodist.py for 
	dissagregation:

Förster, K., Hanzer, F., Winter, B., Marke, T., and Strasser, U.: An open-source 
MEteoroLOgical observation time series DISaggregation Tool 
(MELODIST v0.1.1), Geosci. Model Dev., 9, 2315-2333, doi:10.5194/gmd-9-2315-2016
, 2016.

https://github.com/kristianfoerster/melodist

"""


import matplotlib.pyplot as plt
import melodist
import scipy.stats

import numpy as np
import pandas as pd

import glob
import matplotlib
import xarray as xr
import calendar3 as cal3
#===============================================================================
# 
# MODULE: Subdaily dissagregaration
#
#===============================================================================




#===============================================================================
# Settings
#===============================================================================

# Do we want to run all methods and plot to look at performance?
test_methods=False 
geotop_out=False
fsm_out=True

# standard calenders to interp non-standard cals to in hist and clim period
nc_standard_clim='/home/joel/sim/qmap/test/pyout/aresult/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
nc_standard_hist='/home/joel/sim/qmap/test/pyout/aresult/ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc'
# station attributes
longitude = 9.80988
latitude = 46.82786
timezone = 1
slope  = 0

calibration_period = slice('2000-01-01', '2015-12-31')
validation_period = slice('2016-01-01', '2016-12-31')
plot_period = slice('2016-09-03', '2030-10-13')


#===============================================================================
# model downloads accounting
#===============================================================================
nchistfiles = glob.glob('/home/joel/sim/qmap/test/pyout/'+ "*historical*")
ncrcp26files = glob.glob('/home/joel/sim/qmap/test/pyout/'+ "*rcp26*")
ncrcp85files = glob.glob('/home/joel/sim/qmap/test/pyout/'+ "*rcp85*")
# file names

ncfiles_vec= nchistfiles,ncrcp26files, ncrcp85files

for ncfiles in ncfiles_vec:
	allfiles = [i.split('pyout/', 1)[1] for i in ncfiles] # c. 9000 models, pars, times
	rootfiles = [i.split('day', 1)[0] for i in allfiles]
	modelPars = list(set(rootfiles)) # c.5000 models, par

	models2 = list(set([i.split('_', 1)[1] for i in rootfiles]))  # c. 60 models
	print(len(models2))
	#26                                                                                                                                                                                                                                                                                                               
	#10
	#25
	#total = 61 models downloaded
	# 44 now available complete eg len(nc_complete) not all models have all parameters
#===============================================================================
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
#===============================================================================
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
	


data_obs_hourly = df_interpol.loc[:,('TA','PINT', 'RH', 'ISWR','ILWR', 'VW')]

data_obs_hourly.rename(columns={'TA': 'temp', 'PINT': 'precip', 'RH': 'hum', 
	'ISWR': 'glob', 'ILWR': 'ilwr', 'VW': 'wind'}, inplace=True)



# unit conversions
data_obs_hourly.precip *=1.#(1/3600.)#mm/hr to to kgm-2s-1 
data_obs_hourly.hum *=100. # 0-1 to 0-100

#data_obs_hourly = df2
#data_obs_hourly = melodist.util.drop_incomplete_days(data_obs_hourly)

"""Aggregates data (hourly to daily values) according to the characteristics
of each variable (e.g., average for temperature, sum for precipitation)"""
data_obs_daily = melodist.util.daily_from_hourly(data_obs_hourly) # 

#======================================================================================
# Methods from melodist
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

#===============================================================================



# check for complete nc files

ncfiles = glob.glob('/home/joel/sim/qmap/test/pyout/aresult/'+ '*_TS_ALL_ll.nc')


nc_complete=[]
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
		print (nc)
		print (df.columns)
		nc_complete.append(nc)
		#break
	except:
		#print "no vars"
		#print df.columns
		next



for nc in nc_complete:
	print(nc)
	#nc=nc_complete[9] no leapS
	#nc=nc_complete[10] # 360 day

	# sthis file has a standard calender and used as target in interpolation (hist or clim periods)
	if  'historical' in nc: 
		ds1 = xr.open_dataset(nc_standard_hist, decode_times=True)
	if  'rcp' in nc: 
		ds1 = xr.open_dataset(nc_standard_clim,decode_times=True)

	ds = xr.open_dataset(nc,decode_times=True)

	datetimeindex = ds.indexes["time"]
	if datetimeindex.is_all_dates is False: 
		dateob = datetimeindex[0]
		cal = dateob.calendar
		print(cal)
		# print(dateob.datetime_compatible)

		#ds_conv = convert_calendar(ds, ds1,'year', 'time')  
		ds_interp = cal3.interp_calendar(ds, ds1, 'time')  

		# now we have to fill na in first and last timestep which results from interpolation
		vars = list(ds_interp.data_vars.keys())  
	
		for var in vars:
			# backfill step 0 
			ds_interp[var][0,:,:] =ds_interp[var][1,:,:] 
			# foward fill step -1
			ds_interp[var][-1,:,:] =ds_interp[var][-2,:,:] 

		# this just checks it worked
		datetimeindex = ds_interp.indexes["time"]
		if datetimeindex.is_all_dates is True: 
			#dateob = datetimeindex[0]
			#cal = dateob.calendar
			print("Succesfully converted "+cal+ "to standard calender")
			#print(dateob.datetime_compatible)

			# rename object for next step
			ds = ds_interp
		else:
			"Something went wrong cal3 conversionfailed"

	else:
		print("Standard calender exists")




	# this bit depend on the model and calender i think 
	# datetimeindex = ds.indexes["time"]
	# if datetimeindex.is_all_dates is False: 
	# 	dateob = datetimeindex[0]
	# 	cal = dateob.calendar
	# 	print(cal)
	# 	print(dateob.datetime_compatible)


	# if datetimeindex.is_all_dates is False: # this means we have 365 or 360 day calender


	# 	<p>	# somehow test for cal typ here
	# 	ds.variables['time'] 

	# 	# 'noleapcal' strategy is to just interpolate leap days 29 Feb
	# 	if cal == 'noleap':
	# 		print 'yes'
	# 		start = ds.variables['time'][0].values 
	# 		end = ds.variables['time'][-1].values

	# 		# convert cftime index to datetime object
	# 		# https://stackoverflow.com/questions/55786995/converting-cftime-datetimejulian-to-datetime
	# 		ds_ccaldatetime = ds.indexes['time'].to_datetimeindex()  
	# 		ds_realdatetime =pd.date_range(ds_ccaldatetime[0],ds_ccaldatetime[-1], freq='D')
	# 		ds.reindex(ds_realdatetime, fill_value=0)

	#  	# here cal is 30 days in each year
	# 	if cal == '360_day':
	# 		print 'yes'
	# 		print 'yes'
	# 		start = ds.variables['time'][0].values 
	# 		end = ds.variables['time'][-1].values

	# 		# convert cftime index to datetime object
	# 		# https://stackoverflow.com/questions/55786995/converting-cftime-datetimejulian-to-datetime
	# 		ds_ccaldatetime = ds.indexes['time'].to_datetimeindex()  
	# 		ds_realdatetime =pd.date_range(ds_ccaldatetime[0],ds_ccaldatetime[-1], freq='D')
	# 		ds.reindex(ds_realdatetime, fill_value=0)



	# if datetimeindex.is_all_dates is False:
	# 	try:
	# 		#http://xarray.pydata.org/en/stable/weather-climate.html#non-standard-calendars-and-dates-outside-the-timestamp-valid-range
	# 		datetimeindex = ds.indexes["time"].to_datetimeindex() 
	# 	except:
	# 		print("Non standard time index 360 day cal, cant create datetime index, skipping")
	# 		continue # cant handle 360day claender


			
	df = ds.sel(lon=longitude, lat=latitude, method='nearest').to_dataframe()
	df2 = df.loc[:, ('tas','tasmin' , 'tasmax','pr', 'hurs', 'rsds', 'rlds', 'ps')]
	#compute wind speed and wind direction
	df2.loc[:,'ws'] = np.sqrt(df.uas**2+df.vas**2)
	df2.loc[:,'wd'] = (180 / np.pi) * np.arctan(df.uas/df.vas) + np.where(df.vas>0,180,np.where(df.uas>0,360,0))

	# unit conversions
	df2.loc[:,'pr'] *=3600*24 # kgm-2s-1 (same as mm /s)-> mm/day 
	df2.rename(columns={'tas': 'temp','tasmax': 'tmax','tasmin': 'tmin', 'pr': 'precip', 'hurs': 'hum', 'rsds': 'glob', 'ws': 'wind'}, inplace=True)


	#df2["datetime"] = datetimeindex

	# set index
	df3 =df2.set_index(datetimeindex) 



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
	#stats.calc_precipitation_stats()

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
	
	# FIll  last values that go to NA in interp
	# foward fill step -1
	station.data_disagg['DW'][-1] =station.data_disagg['DW'][-2] 


	# make air pressure equal each hour
	station.data_disagg['P'] = df4.ps.resample('H').mean().interpolate(method='pad') 
	station.data_disagg['P'][-1] =station.data_disagg['P'][-2] 


	# disagg lwin with T

	# compute a daily scaling factor and apply to each hour
	scalingFact = df4.rlds/df4.temp
	scalingFact_day = scalingFact.resample('H').mean()
	scale_fact_day_interp = scalingFact_day.interpolate(method='pad')   
	station.data_disagg['ILWR']  = station.data_disagg.temp * scale_fact_day_interp
	station.data_disagg['ILWR'][-1] =station.data_disagg['ILWR'][-2] 
	#plot(data_obs_hourly.ilwr,station.data_disagg.ILWR)


	# write out netcdf here

	if geotop_out is True:
		outname=nc.split('/')[-1] 
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
		df_gtop.to_csv(path_or_buf="/home/joel/Desktop/"+outname+".txt",na_rep=-999,float_format='%.3f')

	if fsm_out is True:
		outname=nc.split('/')[-1]  
		dates=station.data_disagg.temp.index  
		df_fsm= pd.DataFrame({	

						
						# "Sf":snowTot/(60.*60.), # prate in mm/hr to kgm2/s
						# "Rf":rainTot/(60.*60.), # prate in mm/hr to kgm2/s
						"TA":station.data_disagg.temp, 
						"RH":station.data_disagg.hum,#*0.01, #meteoio 0-1
						"VW":station.data_disagg.wind,
						"DW":station.data_disagg.DW, 
						"P":station.data_disagg.P,
						"ISWR":station.data_disagg.glob, 
						"ILWR":station.data_disagg.ILWR, 
						"PINT":station.data_disagg.precip # prate mm/hr

						})
		df_fsm.index.name="Date"
		dffsm = df_fsm.ffill() # ensure last day of nans (DW,ILWR,P) is filled with last known values (same as other variables, they are just done at the daily level already)
		dffsm.to_csv(path_or_buf='/home/joel/sim/qmap/topoclim/'+outname+'_HOURLY.txt' ,na_rep=-999,float_format='%.8f', header=True, sep=',', 
			columns=['TA', 'RH', 'VW', 'DW', 'P', 'ISWR', 'ILWR','PINT'])
		print("WRITTEN: "+ outname)


		# import pycat

		# from pycat.io import Dataset
		# from pycat.esd import QuantileMapping

		# obs = Dataset('sample-data', 'observation.nc')
		# mod = Dataset('sample-data', 'model*.nc')
		# sce = Dataset('sample-data', 'scenario*.nc')

		# qm = QuantileMapping(obs, mod, sce)


		# qm.correct()
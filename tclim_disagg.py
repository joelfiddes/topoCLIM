import melodist
import pandas as pd
import glob
import numpy as np
cal_period = slice('2000-01-01', '2015-12-31')

def main(daily_cordex,hourly_obs, mylon, mylat,mytz, slope):

	print(mytz)
	df_obs= pd.read_csv(hourly_obs, index_col=0, parse_dates=True)
	data_obs_hourly = df_obs.loc[:,('TA', 'RH', 'ISWR','ILWR', 'VW')]
	data_obs_hourly.rename(columns={'TA': 'temp',  'RH': 'hum', 
		'ISWR': 'glob', 'ILWR': 'ilwr', 'VW': 'wind'}, inplace=True)

	# unit conversions
	data_obs_hourly['precip'] =df_obs.Rf + df_obs.Sf#mm/hr to to kgm-2s-1 
	#data_obs_hourly.hum *=100. # 0-1 to 0-100


	"""Aggregates data (hourly to daily values) according to the characteristics
	of each variable (e.g., average for temperature, sum for precipitation)"""
	data_obs_daily = melodist.util.daily_from_hourly(data_obs_hourly) # 

	# cordex qmap file
	df_cordex= pd.read_csv(daily_cordex, index_col=0, parse_dates=True)


	#compute wind speed and wind direction
	# df2.loc[:,'ws'] = np.sqrt(df.uas**2+df.vas**2)
	# df2.loc[:,'wd'] = (180 / np.pi) * np.arctan(df.uas/df.vas) + np.where(df.vas>0,180,np.where(df.uas>0,360,0))

	# unit conversions
	df_cordex.loc[:,'PINT'] *=3600*24 # kgm-2s-1 (same as mm /s)-> mm/day 
	df_cordex.rename(columns={'TA': 'temp','TAMAX': 'tmax','TAMIN': 'tmin', 'PINT': 'precip', 'RH': 'hum', 'ISWR': 'glob', 'VW': 'wind'}, inplace=True)
	#df3 =df2.set_index(datetimeindex) 
	#df_cordex = df3.resample('D').mean()
	
	# filter tmin /tmax error in #'/home/joel/sim/qmap/test/pyout/aresult/MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
	# vals of 500 / 20
	df_cordex.tmin[df_cordex.tmin<220]=np.nan

	df_cordex.tmax[df_cordex.tmax<220]=np.nan
	df_cordex.temp[df_cordex.temp<220]=np.nan
	df_cordex.temp[df_cordex.temp>330]=np.nan
	# for case of 360 calender that ends on dec 26 (rare) for example:
	#'/home/joel/sim/qmap/test/pyout/aresult/MOHC-HadGEM2-ES_historical_r1i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
	# we run ffill() for last 5 days
	# has no effect if contains data
	df_cordex=df_cordex.ffill()

	station = melodist.Station(lon=float(mylon), lat=float(mylat), timezone=int(mytz) , data_daily=df_cordex)
	station.statistics = melodist.StationStatistics(data_obs_hourly.loc[cal_period])

	stats = station.statistics
	stats.calc_wind_stats()
	stats.calc_humidity_stats()
	stats.calc_temperature_stats()
	stats.calc_radiation_stats()

	#final dissagregation with selected methods
	station.disaggregate_temperature(method='sine_mean')
	station.disaggregate_humidity(method='equal')
	station.disaggregate_radiation(method='pot_rad')
	station.disaggregate_wind(method='random')
	station.disaggregate_precipitation(method='equal')

	# make wind direction equal each hour
	#station.data_disagg['DW'] = df_cordex.DW.resample('H').mean().interpolate(method='pad') 
	
	# FIll  last values that go to NA in interp
	# foward fill step -1
	#station.data_disagg['DW'][-1] =station.data_disagg['DW'][-2] 

	# make air pressure equal each hour
	station.data_disagg['P'] = df_cordex.P.resample('H').mean().interpolate(method='pad') 
	#station.data_disagg['P'][-1] =station.data_disagg['P'][-2] 

	# disagg lwin with T
	# compute a daily scaling factor and apply to each hour
	scalingFact = df_cordex.ILWR/df_cordex.temp
	scalingFact_day = scalingFact.resample('H').mean()
	scale_fact_day_interp = scalingFact_day.interpolate(method='pad')   
	station.data_disagg['ILWR']  = station.data_disagg.temp * scale_fact_day_interp
#	station.data_disagg['ILWR'][-1] =station.data_disagg['ILWR'][-2] 
	#plot(data_obs_hourly.ilwr,station.data_disagg.ILWR)

	# partition prate to rain snow (mm/hr)	
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

	# # linearly reduce snow to zero in steep slopes
	#if steepSnowReduce=="TRUE": # make this an option if need that in future

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



	#Filter TA for absolute 0 vals - still need this?
	station.data_disagg.temp[station.data_disagg.temp<220]=np.nan
	station.data_disagg.temp = station.data_disagg.temp.ffill()
	station.data_disagg.P[station.data_disagg.P<10000]=np.nan
	station.data_disagg.P = station.data_disagg.P.ffill()
	station.data_disagg.ILWR[station.data_disagg.ILWR<100]=np.nan
	station.data_disagg.ILWR = station.data_disagg.ILWR.ffill()
	station.data_disagg.hum[station.data_disagg.hum<5]=np.nan
	station.data_disagg.hum = station.data_disagg.hum.ffill()

	outname=daily_cordex.split('.txt')[0]+  '_F.txt' 
	dates=station.data_disagg.temp.index 

	df_fsm= pd.DataFrame({	

	 				"year": dates.year, 
	 				"month": dates.month, 
					"day": dates.day, 
					"hour": dates.hour,
					"ISWR":station.data_disagg.glob, 
					"ILWR":station.data_disagg.ILWR, 
					"Sf":snowTot/3600, # prate in mm/hr to kgm2/s
					"Rf":rainTot/3600, # prate in mm/hr to kgm2/s
					"TA":station.data_disagg.temp, 
					"RH":station.data_disagg.hum,#*0.01, #meteoio 0-1
					"VW":station.data_disagg.wind,
					"P":station.data_disagg.P,
					
					})
	df_fsm.index.name="Date"
	dffsm = df_fsm.ffill() # ensure last day of nans (DW,ILWR,P) is filled with last known values (same as other variables, they are just done at the daily level already)
	
	#dffsm_3h = dffsm.resample("3H").mean() #dont need know we do inline processing
	dffsm.to_csv(path_or_buf=outname ,na_rep=-999,float_format='%.4f', header=False, sep='\t', index=False, 
		columns=['year','month','day', 'hour', 'ISWR', 'ILWR', 'Sf', 'Rf', 'TA', 'RH', 'VW', 'P'])
	print("WRITTEN: "+ outname)




if __name__ == '__main__':
	import sys
	daily_cordex[1]
	hourly_obs = sys.argv[2]
	mylon = sys.argv[3]
	mylat = sys.argv[4]
	mytz = sys.argv[5]
	main(daily_cordex,hourly_obs,  mylon, mylat,mytz)			
import melodist
import pandas as pd
import glob
import numpy as np
cal_period = slice('2000-01-01', '2015-12-31')

def main(daily_cordex,hourly_obs, mylon, mylat,mytz):

	print(mytz)
	df_obs= pd.read_csv(hourly_obs, index_col=0, parse_dates=True)
	data_obs_hourly = df_obs.loc[:,('TA','PINT', 'RH', 'ISWR','ILWR', 'VW')]
	data_obs_hourly.rename(columns={'TA': 'temp', 'PINT': 'precip', 'RH': 'hum', 
		'ISWR': 'glob', 'ILWR': 'ilwr', 'VW': 'wind'}, inplace=True)

	# unit conversions
	data_obs_hourly.precip *=(1/3600.)#mm/hr to to kgm-2s-1 
	data_obs_hourly.hum *=100. # 0-1 to 0-100

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
	station.disaggregate_radiation(method='pot_rad_via_bc')
	station.disaggregate_wind(method='random')
	station.disaggregate_precipitation(method='equal')

	# make wind direction equal each hour
	station.data_disagg['DW'] = df_cordex.DW.resample('H').mean().interpolate(method='pad') 
	
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

	outname=daily_cordex.split('.txt')[0]+  '_H.txt' 
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
	dffsm.to_csv(path_or_buf=outname ,na_rep=-999,float_format='%.8f', header=True, sep=',', 
		columns=['TA', 'RH', 'VW', 'DW', 'P', 'ISWR', 'ILWR','PINT'])
	print("WRITTEN: "+ outname)




if __name__ == '__main__':
	import sys
	daily_cordex[1]
	hourly_obs = sys.argv[2]
	mylon = sys.argv[3]
	mylat = sys.argv[4]
	mytz = sys.argv[5]
	main(daily_cordex,hourly_obs,  mylon, mylat,mytz)			
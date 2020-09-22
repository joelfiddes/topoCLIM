# a more evolved way of tclim' ing

"""
python3

This module takes post processed daily CORDEX downloads from esgf_post.py and 
produces point timeseries with standard calenders

Example:


Vars:


Details:


"""
from configobj import ConfigObj
import glob
import xarray as xr
import calendar3 as cal3

wd = "/home/joel/sim/qmap"
raw_dir= wd+'/raw_cordex/'
nc_standard_clim=raw_dir+'/aresult/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
nc_standard_hist=raw_dir+'/aresult/ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc'


def howmanyFiles(wd)
	config = ConfigObj(wd+"/config.ini")


	#===============================================================================
	# model downloads accounting
	#===============================================================================
	nchistfiles = glob.glob(raw_dir+ "*historical*")
	ncrcp26files = glob.glob(raw_dir+ "*rcp26*")
	ncrcp85files = glob.glob(raw_dir+ "*rcp85*")
	# file names

	ncfiles_vec= nchistfiles,ncrcp26files, ncrcp85files

	for ncfiles in ncfiles_vec:
		allfiles = [i.split('raw_cordex/', 1)[1] for i in ncfiles] # c. 9000 models, pars, times
		rootfiles = [i.split('day', 1)[0] for i in allfiles]
		modelPars = list(set(rootfiles)) # c.5000 models, par

		models2 = list(set([i.split('_', 1)[1] for i in rootfiles]))  # c. 60 models
		print(len(models2))
		#26                                                                                                                                                                                                                                                                                                               


	# check for complete nc files
def completeFiles(raw_dir):
	ncfiles = glob.glob(raw_dir+'/aresult/'+ '*_TS_ALL_ll.nc')


	nc_complete=[]
	for nc in ncfiles:

		ds = xr.open_dataset(nc,decode_times=True)

		if  len(ds.variables.keys())  <15:
			print ("not enough vars "+nc)
			next   
		# try:

		# 	df2 = df[['tas','tasmin' , 'tasmax','pr', 'hurs', 'rsds', 'rlds', 'ps']]#, 'vas']]  
		# 	#print 'ok'
		# 	print (nc)
		# 	print (df.columns)
		# 	nc_complete.append(nc)
			#break
		#except:
		else:
			print (nc)
			nc_complete.append(nc)

	return(nc_complete)


def calendarNinja(nc_complete, nc_standard_hist,nc_standard_clim):
	for nc in nc_complete:
		print(nc)
		outname = nc.split(".nc")[0]+"_SCAL.nc" 
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
				ds.to_netcdf(outname)  
			else:
				"Something went wrong cal3 conversionfailed"

		else:
			print("Standard calender exists " + nc)
			os.rename(nc, outname)


nc_complete = completeFiles(raw_dir)  
calendarNinja(nc_complete, nc_standard_hist,nc_standard_clim)




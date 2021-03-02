require(qmap)
require(ncdf4)
#require(caTools)
##require(pracma)
#require(hydroGOF)
#require(wesanderson)
# Major rework: 12 Oct 2020
# Args:
	# indir="/home/joel/sim/qmap/topoclim_test" 

# Example:
	# Rscript', 'qmap_hour_plots_daily.R', '/home/joel/sim/qmap/tclim_points/', 'tscale_1_1D', '/home/joel/sim/qmap/tscale_points//out/tscale_1_1D.csv', '9.809', '46.83', '/home/joel/sim/qmap/raw_cordex/aresult/'

args = commandArgs(trailingOnly=TRUE)
wd = args[1]
sample = args[2]
daily_obs = args[3]
mylon = as.numeric(args[4])
mylat = as.numeric(args[5])
CORDEXPATH = args[6]

# wd='/home/joel/sim/qmap/tclim_points/'
# sample = "tscale_1_1D"
# daily_obs ='/home/joel/sim/qmap/tscale_points//out/tscale_1_1D.csv'
# mylon = 9.809
# mylat = 46.83
# CORDEXPATH='/home/joel/sim/qmap/raw_cordex/aresult/'

# linearly resample 3h era5 to 1h and cpture path as variable
#daily_obs = system2(command = "python", args = paste( "resample_timeseries.py" , path_inp), stdout=TRUE)

		

# Setup ========================================================================
indir=paste0(wd, '/s',sample)
dir.create(indir)

outdir=paste0(indir,"/aqmap_results/")
dir.create(outdir)

# dissaggregated calender corrected hourly input files 
files = list.files(path=CORDEXPATH, pattern="SCAL.nc", recursive=F, full.name=T)
hist_files = files[ grep('historical',files)]
rcp26_files = files[ grep('rcp26',files)]
rcp85_files = files[ grep('rcp85',files)]

if ( length(hist_files)==0){print("ERROR! No cordex files found!")}
# different models have different time spans so need to define common periods where all models have data

# common historical period of models
startDateHist = "1979-01-01"
endDateHist = "2005-12-31"

# common rcp periods of models - perhaps these can/should be 1 set, in practice is the same
startDateRcp26 = "2006-01-02"
endDateRcp26 = "2099-12-31"
startDateRcp85 = "2006-01-02"
endDateRcp85 = "2099-12-31"

# the period to do quantile mapping over (obs and sim must overlap!)
startDateQmap = "1980-01-01"
endDateQmap = "2004-12-31"





# Obs Era5 data downscaled by toposcale, can be resampled

obs=read.csv(daily_obs)
myvars=c('tas', 'tasmin', 'tasmax','pr', 'ps', 'hurs', 'rsds','rlds', 'uas')

OItas = which(names(obs)=='TA')
OItasmin = which(names(obs)=='TAMIN')
OItasmax = which(names(obs)=='TAMAX')
OISf = which(names(obs)=='Sf')
OIRf = which(names(obs)=='Rf')
OIp = which(names(obs)=='P')
OIrh = which(names(obs)=='RH')
OIiswr = which(names(obs)=='ISWR')
OIilwr = which(names(obs)=='ILWR')
OIvw = which(names(obs)=='VW')


# names(obs)
#  [1] "datetime" "ISWR"     "ILWR"     "Sf"       "Rf"       "TA"      
#  [7] "RH"       "VW"       "P"        "TAMAX"    "TAMIN" 

for (var in myvars){
print(var)



	if(var=="tas"){obsindex <-OItas; convFact <-1} # K
	if(var=="tasmin"){obsindex <-OItasmax; convFact <-1} # K
	if(var=="tasmax"){obsindex <-OItasmin; convFact <-1} # K
		if(var=="pr"){obsindex <-OISf; convFact <-(1)} # both FSm and cordex are kgm-2s-1 
			if(var=="ps"){obsindex <-OIp; convFact <-1}# Pa
				if(var=="hurs"){obsindex <-OIrh; convFact <-1} # % 0-100
					if(var=="rsds"){obsindex <-OIiswr; convFact <-1}	# Wm-2
						if(var=="rlds"){obsindex <-OIilwr; convFact <-1}# Wm-2
							if(var=="uas"){obsindex <-OIvw; convFact <-1}# ms-1


	# read and convert obs unit
	obs_var=obs[,obsindex]*convFact	
	
	# if pr need to add obs Sf amnd Rf components
	if (var=='pr'){
	obs_var=(obs[,OIRf]+obs[,OISf])*convFact	
	}

	# aggregate obs to daily resolution
	obs_datetime= strptime(obs$datetime,format="%Y-%m-%d")

	# cut obs to qmap period 
	start_obs = which(obs_datetime==startDateQmap)
	end_obs=which(obs_datetime==endDateQmap)
	obs_qmap_period= obs_var[start_obs: end_obs]

	# cut obs to cordex historical common  period 
	start_obs =1 # simply first ob as CORDEX hist period stats before 1980-01-01 otherwise use "end_obs=which(obs_datetime==startDateHist)"
	end_obs=which(obs_datetime==endDateHist)
	obs_cp_period= obs_var[start_obs: end_obs]
	obs_cp_dates = obs_datetime[start_obs: end_obs]


	
	# Cordex historical "1970-01-01 12" UTC" to 2005-12-30 12"======================
	hist_qmap_list <- list()
	hist_nqmap <- list()
	hist_dates<-list()
	hist_qmap_season_list <- list()
	modelChain_vec=c()

	for (hist_file in hist_files){
	# check if done
	


		print(hist_file)
		modelNameBase = unlist(strsplit(hist_file,'/'))[length(unlist(strsplit(hist_file,'/')))]
		GCM = unlist(strsplit(modelNameBase,'_historical'))[1]
		RCM = unlist(strsplit(modelNameBase,'_'))[4]
		modelChain=paste0(GCM,'_',RCM)
		
		
		# get data READ netcdf

		nc=nc_open(hist_file)
		tas_allgrids =ncvar_get(nc, var)


		lon =ncvar_get(nc, 'lon')
		lat =ncvar_get(nc, 'lat')
		time =ncvar_get(nc, 'time')
		myx =which.min(abs(lon - mylon))
		myy = which.min(abs(lat - mylat))

		# extract timeseries corresponding to qmap (all)
		tas = tas_allgrids[myx,myy,]

		# wspeed compute from u and v here
		if(var=="uas"){
			vas_allgrids =ncvar_get(nc, 'vas')
			vas = tas_allgrids[myx,myy,]
			tas = sqrt(tas^2 + vas^2)
		  	}





		# define time
		time = ncvar_get( nc,'time')
		time2=time
		z <- time2*24*60*60 #convert days to seconds
		origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
		origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section
		datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
		cordex_dates  = as.Date(datesPl)


		if(any(is.na(tas)==TRUE)){print(paste0("no data found in variable, skipping this file", hist_file));next}
		
		# if all present and correct can now add to chain
		modelChain_vec = c(modelChain_vec,modelChain) # available hist data

		# define periods

		start_cp_cord = which(strptime(cordex_dates, format="%Y-%m-%d")==startDateHist)
		end_cp_cord=which( strptime(cordex_dates, format="%Y-%m-%d")==endDateHist  )
		start_qmap_cord = which(strptime(cordex_dates, format="%Y-%m-%d")==startDateQmap)
		end_qmap_cord=which( strptime(cordex_dates, format="%Y-%m-%d")==endDateQmap  )

		
		# extract data of common period historical timeseries 1979-2005 
		hist_cp=tas[start_cp_cord:end_cp_cord]

		# extract dates common period historical timeseries 1979-2005 
		cordex_dates_cp = cordex_dates[start_cp_cord:end_cp_cord]

		# extract dat for qmap period
		hist_qmap_period=tas[start_qmap_cord:end_qmap_cord]

		# get gmap pars
		pars = fitQmap(obs_qmap_period,hist_qmap_period, method = "QUANT")

		# do qmap
		hist_qmapped = doQmap(hist_cp,pars)

		# add to list
		hist_qmap_list[[modelChain]] <- list(hist_qmapped)
		hist_nqmap[[modelChain]] <- list(hist_cp)

		# do seasonal qmap
		month_cp = substring(cordex_dates_cp,6,7)
		summer_index=which( month_cp=="06"| month_cp=="07"| month_cp=="08")
		autumn_index=which(month_cp=="09"| month_cp=="10"| month_cp=="11")
		winter_index=which(month_cp=="12"| month_cp=="01"| month_cp=="02")
		spring_index=which(month_cp=="03"| month_cp=="04"| month_cp=="05")
		
		pars_summer = fitQmap(obs_qmap_period[summer_index],hist_qmap_period[summer_index], method = "QUANT")
		pars_autumn = fitQmap(obs_qmap_period[autumn_index],hist_qmap_period[autumn_index], method = "QUANT")
		pars_winter = fitQmap(obs_qmap_period[winter_index],hist_qmap_period[winter_index], method = "QUANT")
		pars_spring = fitQmap(obs_qmap_period[spring_index],hist_qmap_period[spring_index], method = "QUANT")

		hist_summer = doQmap(hist_cp[summer_index],pars_summer)
		hist_autumn = doQmap(hist_cp[autumn_index],pars_autumn)
		hist_winter= doQmap(hist_cp[winter_index],pars_winter)
		hist_spring = doQmap(hist_cp[spring_index],pars_spring)

		hist_season<-1:length(hist_cp)
		hist_season[summer_index]<-hist_summer
		hist_season[autumn_index]<-hist_autumn
		hist_season[winter_index]<-hist_winter
		hist_season[spring_index]<-hist_spring

		hist_qmap_season_list[[modelChain]] <- list(hist_season)
		hist_qmap_list[[modelChain]][['pars']]<-pars 
		hist_qmap_list[[modelChain]][['pars_summer']]<-pars_summer
		hist_qmap_list[[modelChain]][['pars_autumn']]<-pars_autumn 
		hist_qmap_list[[modelChain]][['pars_winter']]<-pars_winter 
		hist_qmap_list[[modelChain]][['pars_spring']]<-pars_spring 
		}

	hist_qmap_list_NOPARS = lapply(hist_qmap_list, `[[`, 1)
	mydata=hist_qmap_list_NOPARS
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-cordex_dates_cp 
	save(df,file = paste0(outdir,"hist_qmap_",var))

	hist_qmap_list_NOPARS = lapply(hist_qmap_season_list, `[[`, 1)
	mydata=hist_qmap_list_NOPARS
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-cordex_dates_cp 
	save(df,file = paste0(outdir,"hist_season_",var))

	mydata=hist_nqmap
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-cordex_dates_cp 
	save(df,file = paste0(outdir,"hist_nmap_",var))

	# Cordex rcp26 =======================================================================

	rcp26_qmap_list <- list()
	rcp26_nqmap <- list()
	rcp26_qmap_season_list <- list()

	for (rcp26_file in rcp26_files){
		# check if done
	
		print(rcp26_file)
		modelNameBase = unlist(strsplit(rcp26_file,'/'))[length(unlist(strsplit(rcp26_file,'/')))]
		GCM = unlist(strsplit(modelNameBase,'_rcp26_'))[1]
		RCM = unlist(strsplit(modelNameBase,'_'))[4]
		modelChain=paste0(GCM,'_',RCM)

		# skip to next model if does not have hist data
		if (!modelChain %in% modelChain_vec ){print(paste0("No historical data for: ", modelChain));next}



		# get data READ netcdf
		nc=nc_open(rcp26_file)
		tas_allgrids =ncvar_get(nc, var) 
		lon =ncvar_get(nc, 'lon')
		lat =ncvar_get(nc, 'lat')
		time =ncvar_get(nc, 'time')
		myx =which.min(abs(lon - mylon))
		myy = which.min(abs(lat - mylat))

		# extract timeseries corresponding to qmap (all)
		tas = tas_allgrids[myx,myy,]

		# wspeed compute from u and v here
		if(var=="uas"){
			vas_allgrids =ncvar_get(nc, 'vas')
			vas = tas_allgrids[myx,myy,]
			tas = sqrt(tas^2 + vas^2)
		  	}

		# define time
		time = ncvar_get( nc,'time')
		time2=time
		z <- time2*24*60*60 #convert days to seconds
		origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
		origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section
		datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
		greg_cal_rcp  = as.Date(datesPl)



		# extract variable

		if(any(is.na(tas)==TRUE)){print(paste0("no data found in variable, skipping this file", rcp26_file));next}



		# do qmap
		rcp26_qmapped = doQmap(tas, hist_qmap_list[[modelChain]][['pars']] )

		# trim to common length
		si = which(greg_cal_rcp==startDateRcp26)
		ei = which(greg_cal_rcp==endDateRcp26)

		# assign to list
		rcp26_qmap_list[[modelChain]] <- list(rcp26_qmapped[si:ei])
		rcp26_nqmap[[modelChain]] <- list(tas[si:ei])

		# create seasonal indexes
		month_cp = substring(greg_cal_rcp,6,7)
		summer_index=which( month_cp=="06"| month_cp=="07"| month_cp=="08")
		autumn_index=which(month_cp=="09"| month_cp=="10"| month_cp=="11")
		winter_index=which(month_cp=="12"| month_cp=="01"| month_cp=="02")
		spring_index=which(month_cp=="03"| month_cp=="04"| month_cp=="05")

		# do seasonal qmap
		rcp26_summer = doQmap(tas[summer_index], hist_qmap_list[[modelChain]][['pars_summer']] )
		rcp26_autumn = doQmap(tas[autumn_index], hist_qmap_list[[modelChain]][['pars_autumn']] )
		rcp26_winter= doQmap(tas[winter_index], hist_qmap_list[[modelChain]][['pars_winter']] )
		rcp26_spring = doQmap(tas[spring_index], hist_qmap_list[[modelChain]][['pars_winter']])

		rcp26_season<-1:length(rcp26_qmapped)
		rcp26_season[summer_index]<-rcp26_summer
		rcp26_season[autumn_index]<-rcp26_autumn
		rcp26_season[winter_index]<-rcp26_winter
		rcp26_season[spring_index]<-rcp26_spring

		rcp26_qmap_season_list[[modelChain]] <- list(rcp26_season[si:ei])
		}
		

	mydata=rcp26_qmap_list
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<- greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp26_qmap_",var))

	mydata=rcp26_qmap_season_list
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp26_season_",var))

	mydata=rcp26_nqmap
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp26_nmap_",var))
	# Cordex rcp85 2006-2100=======================================================================

	rcp85_qmap_list <- list()
	rcp85_nqmap <- list()
	rcp85_qmap_season_list <- list()
	for (rcp85_file in rcp85_files){
		# check if done
	
		print(rcp85_file)

		modelNameBase = unlist(strsplit(rcp85_file,'/'))[length(unlist(strsplit(rcp85_file,'/')))]
		GCM = unlist(strsplit(modelNameBase,'_rcp85_'))[1]
		RCM = unlist(strsplit(modelNameBase,'_'))[4]
		modelChain=paste0(GCM,'_',RCM)

		# kip to next model if does not have hist data
		if (!modelChain %in% modelChain_vec ){print(paste0("No historical data for: ", modelChain));next}

		# get data READ netcdf
		nc=nc_open(rcp85_file)
		tas_allgrids =ncvar_get(nc, var) 
		lon =ncvar_get(nc, 'lon')
		lat =ncvar_get(nc, 'lat')
		time =ncvar_get(nc, 'time')
		myx =which.min(abs(lon - mylon))
		myy = which.min(abs(lat - mylat))

		# extract timeseries corresponding to qmap (all)
		tas = tas_allgrids[myx,myy,]

		# wspeed compute from u and v here
		if(var=="uas"){
			vas_allgrids =ncvar_get(nc, 'vas')
			vas = tas_allgrids[myx,myy,]
			tas = sqrt(tas^2 + vas^2)
		  	}

		# define time
		time = ncvar_get( nc,'time')
		time2=time
		z <- time2*24*60*60 #convert days to seconds
		origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
		origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section
		datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
		greg_cal_rcp  = as.Date(datesPl)

		if(any(is.na(tas)==TRUE)){print(paste0("no data found in variable, skipping this file", rcp85_file));next}

	

		# do qmap
		rcp85_qmapped = doQmap(tas, hist_qmap_list[[modelChain]][['pars']] )

		# trim to common length
		si = which(greg_cal_rcp==startDateRcp26)
		ei = which(greg_cal_rcp==endDateRcp26)

		# assign to list
		rcp85_qmap_list[[modelChain]] <- list(rcp85_qmapped[si:ei])
		rcp85_nqmap[[modelChain]] <- list(tas[si:ei])
	
		# create seasonal indexes
		month_cp = substring(greg_cal_rcp,6,7)
		summer_index=which( month_cp=="06"| month_cp=="07"| month_cp=="08")
		autumn_index=which(month_cp=="09"| month_cp=="10"| month_cp=="11")
		winter_index=which(month_cp=="12"| month_cp=="01"| month_cp=="02")
		spring_index=which(month_cp=="03"| month_cp=="04"| month_cp=="05")

		# do seasonal qmap
		rcp85_summer = doQmap(tas[summer_index], hist_qmap_list[[modelChain]][['pars_summer']] )
		rcp85_autumn = doQmap(tas[autumn_index], hist_qmap_list[[modelChain]][['pars_autumn']] )
		rcp85_winter= doQmap(tas[winter_index], hist_qmap_list[[modelChain]][['pars_winter']] )
		rcp85_spring = doQmap(tas[spring_index], hist_qmap_list[[modelChain]][['pars_winter']])

		rcp85_season<-1:length(rcp85_qmapped)
		rcp85_season[summer_index]<-rcp85_summer
		rcp85_season[autumn_index]<-rcp85_autumn
		rcp85_season[winter_index]<-rcp85_winter
		rcp85_season[spring_index]<-rcp85_spring

		rcp85_qmap_season_list[[modelChain]] <- list(rcp85_season[si:ei])
		}


	mydata=rcp85_qmap_list
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<- greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp85_qmap_",var))

	mydata=rcp85_qmap_season_list
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp85_season_",var))

	mydata=rcp85_nqmap
	modelNames = names(mydata)
	df=as.data.frame(mydata)
	names(df)<-modelNames
	df$Date<-greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp85_nmap_",var))

}



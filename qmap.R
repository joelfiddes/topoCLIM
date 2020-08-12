require(qmap)
require(ncdf4)
require(caTools)
require(PCICt) # converts to 365 day clim model calender
require(pracma)
require(hydroGOF)
require(wesanderson)

# TAS script
# handle different calender
# Setup ========================================================================

# POI (WFJ) - to extract cordex gridbox
mylon = 9.80490
mylat = 46.83066

# set end date of climate timeseires to cut to as different model provide differnet length timeseries
ENDDATE_CLIM="2099-12-31"

# files
indir="/home/joel/sim/qmap/test/pyout/aresult"
outdir=paste0(indir,"/aqmap_results/")
files = list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)
dir.create(outdir)
hist_files = files[ grep('historical',files)]
rcp26_files = files[ grep('rcp26',files)]
rcp85_files = files[ grep('rcp85',files)]

# define periods
startDateHist = "1979-01-01"
endDateHist = "2005-12-31"
startDateRcp26 = "2006-01-02"
endDateRcp26 = "2099-12-31"
startDateRcp85 = "2006-01-02"
endDateRcp85 = "2099-12-31"

# Cordex data
#hist_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/historical/r1i1p1/REMO2015/v1/day/tas/tas_EUR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"
#rcp26_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/rcp26/r1i1p1/REMO2015/v1/day/tas/tas_EUR-22_MOHC-HadGEM2-ES_rcp26_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"
#rcp85_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/rcp85/r1i1p1/REMO2015/v1/day/tas/tas_EUR-22_MOHC-HadGEM2-ES_rcp85_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"

pdf(paste(outdir,"evalplot_singleModel_season.pdf"),height=20, width=10)
par(mfrow=c(8,2))

# Obs Era5 data downscaled by toposcale, can be resampled
obsfile = "/home/joel/sim/qmap/wfj_long.csv"#
myvars=c('tas', 'pr', 'ps', 'hurs', 'rsds','rlds', 'uas', 'vas')
for (var in myvars){
print(var)
	# Obs 1980-01-01 - 2017-12-31 ==========================================================================
	#read obs variable TA


	if(var=="tas"){obsindex <-11; convFact <-1} # K
		if(var=="pr"){obsindex <-6; convFact <-(1/3600)} # kgm-2s-1
			if(var=="ps"){obsindex <-5; convFact <-1}# Pa
				if(var=="hurs"){obsindex <-8; convFact <-100} # % 0-100
					if(var=="rsds"){obsindex <-4; convFact <-1}	# Wm-2
						if(var=="rlds"){obsindex <-3; convFact <-1}# Wm-2
							if(var=="uas"){obsindex <-2; convFact <-1}# ms-1
								if(var=="vas"){obsindex <-2; convFact <-1}# ms-1

	obs=read.csv(obsfile)
	obs_var=obs[,obsindex]*convFact	# !!HARDCODE!!

	# compute u and v here
	if(var=="uas"){
	ws <- obs$VW
	wd <- obs$DW
	    obs_var = ws * cos(wd)
	  }
	    if(var=="vas"){
	    ws <- obs$VW
		wd <- obs$DW
    	obs_var = ws * sin(wd)
    }
	# aggregate obs to daily resolution
	days=substr(obs$datetime,1,10)
	d =aggregate(obs_var, list(days), 'mean')
	obs_daily=d$x
	obs_daily_date=d$Group.1


	# Cordex historical "1970-01-01 12" UTC" to 2005-12-30 12"======================
	hist_qmap_list <- list()
	hist_nqmap <- list()
	hist_dates<-list()
	hist_qmap_season_list <- list()

	for (hist_file in hist_files){

		#modfile <- hist_files[1]
		# assign stuff
		modfile=hist_file

		# define variable and gridbox
		nc=nc_open(modfile)
		if(!var %in% names(nc$var) ){print("not enough vars, skipping this file");print(length(nc$var));next}
	   
		 		# !!HARDCODE!!
		tas =ncvar_get(nc, var) 

		lon =ncvar_get(nc, 'lon')
		lat =ncvar_get(nc, 'lat')
		myx =which.min(abs(lon - mylon))
		myy = which.min(abs(lat - mylat))

		# extract timeseries corresponding to qmap (all)
		tas_all = tas[myx,myy,]
		if(all(is.na(tas_all)==TRUE)){print("no data found, skipping this file");next}
		if(anyNA(tas_all)==TRUE){print("NAs found, skipping this file");next}
		# define time
		time = ncvar_get( nc,'time')
		time2=time
		z <- time2*24*60*60 #convert days to seconds
		#origin=substr(nc$dim$time$units,13,20)
		origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
		origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section

		# calender dependent section
		cal = nc$dim$time$calendar
		print(cal)
		# standard cal is assumed to equal proleptic_gregorian
		if(cal=='standard'){cal <- "proleptic_gregorian"}
		if(cal=='gregorian'){cal <- "proleptic_gregorian"}

		if(cal=="360_day"){
			cal <- "360_day"
			seconds.per.day <- 86400
			ts.dat.days <- time
			origin.pcict <- as.PCICt(origin, cal)
			ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
			datesPl <- ts.dat.pcict
			}
		if(cal=="365_day"){
			cal <- "365_day"
			seconds.per.day <- 86400
			ts.dat.days <- time
			origin.pcict <- as.PCICt(origin, cal)
			ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
			datesPl <- ts.dat.pcict
			}

		if(cal=="proleptic_gregorian"){
			datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
			}

		# convert 360 to 365/366 calender with interpolation
		#greg_cal = seq(as.Date('1970-01-01'), as.Date('2005-12-31'), "days") # hardcoded
		if(cal=="360_day"){
			greg_cal_hist = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10))+1, "days") # add plus 1 here as cordex always ends on dec 30
			}

		if(cal=="365_day"){
			greg_cal_hist = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") 
			}

		if(cal=="proleptic_gregorian"){
			greg_cal_hist = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") # no +1 need

			# catch exmple where greg cal ends on dec 30: "/home/joel/sim/qmap/test/pyout/aresult/MPI-M-MPI-ESM-LR_historical_r1i1p1_CLMcom-CCLM4-8-17_v1__TS.nc_TS_ALL_ll.nc"
			if(substring(datesPl[length(datesPl)],1,10)=="2005-12-30"){
					greg_cal_hist = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10))+1, "days") # no +1 need
			}}

		# get index for interp1d
		greg_unix = as.numeric(as.POSIXct(greg_cal_hist, format="%Y-%m-%d"))
		cordCal_unix = as.numeric(as.POSIXct( as.Date(substr(datesPl,1,10)), format="%Y-%m-%d"))

		# NAs generated in cordCal_unix as dates are made in the pcict calender conversion that dont exists ie Feb 30
		# we remove these values completely representing 0.48% of values (63/12960)
		vals2rm = which(is.na(cordCal_unix)==T)
		vals2keep = which(!is.na(cordCal_unix)==T)


		overlap_index = which(greg_unix %in% cordCal_unix)
		tas_365_minus1 = interp1(x=cordCal_unix[vals2keep], y=tas_all[vals2keep], xi = greg_unix[1:(length(greg_unix)-1)], method = c("linear"))

		# we miss very last dec 31 as x only goes to dec 30 and inter1 cannot interp beyond range of x (extrapolate). thats why we do this -1: xi = greg_unix[1:(length(greg_unix)-1)]
		# to fix we simply repeat last value
		tas_365 = c(tas_365_minus1, tas_365_minus1[length(tas_365_minus1)])

		#datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
		days_cordex=greg_cal_hist
		years_hist=as.numeric(substr(greg_cal_hist,1,4))
		start_cp_cord = which(days_cordex==days[1])
		end_cp_cord=which(days_cordex==days[length(days)])

		if(length(end_cp_cord)==0){ # if true this mean that obs extend beyond end_cp_cord of the cordex data - this is the case for cordex historical
			end_cp_cord= length(days_cordex) # then simply take last cordex date
			}
		if(length(start_cp_cord)==0){ # if true this mean that obs starts before cordex data - this is the case for cordex rcp data
			start_cp_cord = 1 # then simply take first cordex date
			}
		# extract timeseries corresponding to obs (subset in time/space)
		hist_cp=tas_365[start_cp_cord:end_cp_cord]

		# define common period
		start_cp= greg_cal_hist[start_cp_cord]
		end_cp = greg_cal_hist[end_cp_cord]

		# cut obs to common period cp
		start_obs = which(obs_daily_date==start_cp)
		end_obs=which(obs_daily_date==end_cp)
		obs_daily_cp = obs_daily[start_obs: end_obs]

		# get gmap pars
		pars = fitQmap(obs_daily_cp,hist_cp, method = "QUANT")

		# do qmap
		hist = doQmap(tas_365,pars)

		#annual aggregated timeseries
		hist_year=aggregate(hist, list(years_hist), 'mean')
		tas_365_year=aggregate(tas_365, list(years_hist), 'mean')
		obs_year=aggregate(obs_daily, list(substring(obs_daily_date,1,4)), 'mean')

		# trim to common length
		si = which(greg_cal_hist==startDateHist)
		ei = which(greg_cal_hist==endDateHist)

		# add to list
		hist_qmap_list[[hist_file]] <- list(hist[si:ei])
		hist_nqmap[[hist_file]] <- list(tas_365[si:ei])
		hist_dates[[hist_file]]<- list(greg_cal_hist[si:ei] )
		


		# do seasonal qmap
		month_cp = substring(greg_cal_hist,6,7)
		summer_index=which( month_cp=="06"| month_cp=="07"| month_cp=="08")
		autumn_index=which(month_cp=="09"| month_cp=="10"| month_cp=="11")
		winter_index=which(month_cp=="12"| month_cp=="01"| month_cp=="02")
		spring_index=which(month_cp=="03"| month_cp=="04"| month_cp=="05")
		
		pars_summer = fitQmap(obs_daily_cp[summer_index],hist_cp[summer_index], method = "QUANT")
		pars_autumn = fitQmap(obs_daily_cp[autumn_index],hist_cp[autumn_index], method = "QUANT")
		pars_winter = fitQmap(obs_daily_cp[winter_index],hist_cp[winter_index], method = "QUANT")
		pars_spring = fitQmap(obs_daily_cp[spring_index],hist_cp[spring_index], method = "QUANT")

		
		hist_summer = doQmap(tas_365[summer_index],pars_summer)
		hist_autumn = doQmap(tas_365[autumn_index],pars_autumn)
		hist_winter= doQmap(tas_365[winter_index],pars_winter)
		hist_spring = doQmap(tas_365[spring_index],pars_spring)

		hist_season<-1:length(hist)
		hist_season[summer_index]<-hist_summer
		hist_season[autumn_index]<-hist_autumn
		hist_season[winter_index]<-hist_winter
		hist_season[spring_index]<-hist_spring


		hist_qmap_season_list[[hist_file]] <- list(hist_season[si:ei])
}

		
		mydata=hist_qmap_list
		mydates=hist_dates
		modelNames = sapply(strsplit(names(mydata),'/'), "[", 9)
		df=as.data.frame(mydata)
		names(df)<-modelNames
		df$Date<-as.Date(unlist(mydates[[1]] ))
		save(df,file = paste0(outdir,"hist_",var))

	# Cordex rcp26 =======================================================================

	rcp26_qmap_list <- list()
	rcp26_nqmap <- list()
	rcp26_dates <-list()
	rcp26_qmap_season_list <- list()
	for (rcp26_file in rcp26_files){

		modfile=rcp26_file
		# define variable and gridbox
		nc=nc_open(modfile)
		

		if(!'tas' %in% names(nc$var) ){print("not enough vars, skipping this file");print(length(nc$var));next}

		tas =ncvar_get(nc, 'tas') 		# !!HARDCODE!!
		lon =ncvar_get(nc, 'lon')
		lat =ncvar_get(nc, 'lat')
		myx =which.min(abs(lon - mylon))
		myy = which.min(abs(lat - mylat))

		# extract timeseries corresponding to qmap (all)
		tas_all = tas[myx,myy,]
		if(all(is.na(tas_all)==TRUE)){print("no data found, skipping this file");next}
		if(anyNA(tas_all)==TRUE){print("NAs found, skipping this file");next}

		# define time
		time = ncvar_get( nc,'time')
		time2=time
		z <- time2*24*60*60 #convert days to seconds
		#origin=substr(nc$dim$time$units,13,20)
		origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
		origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section

		# calender dependent section
		cal = nc$dim$time$calendar
		print(cal)
		# standard cal is assumed to equal proleptic_gregorian
		if(cal=='standard'){cal <- "proleptic_gregorian"}
		if(cal=='gregorian'){cal <- "proleptic_gregorian"}

		if(cal=="360_day"){
			cal <- "360_day"
			seconds.per.day <- 86400
			ts.dat.days <- time
			origin.pcict <- as.PCICt(origin, cal)
			ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
			datesPl <- ts.dat.pcict
			}
		if(cal=="365_day"){
			cal <- "365_day"
			seconds.per.day <- 86400
			ts.dat.days <- time
			origin.pcict <- as.PCICt(origin, cal)
			ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
			datesPl <- ts.dat.pcict
			}

		if(cal=="proleptic_gregorian"){
			datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
			}

		# convert 360 to 365/366 calender with interpolation
		#greg_cal = seq(as.Date('1970-01-01'), as.Date('2005-12-31'), "days") # hardcoded
		if(cal=="360_day"){
			greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(ENDDATE_CLIM), "days")
			}

		if(cal=="365_day"){
			greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") 
			}

		if(cal=="proleptic_gregorian"){
			greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(ENDDATE_CLIM), "days") # no +1 need
			}


		# get index for interp1d
		greg_unix = as.numeric(as.POSIXct(greg_cal_rcp, format="%Y-%m-%d"))
		cordCal_unix = as.numeric(as.POSIXct( as.Date(substr(datesPl,1,10)), format="%Y-%m-%d"))

		# NAs generated in cordCal_unix as dates are made in the pcict calender conversion that dont exists ie Feb 30
		# we remove these values completely representing 0.48% of values (63/12960)
		vals2rm = which(is.na(cordCal_unix)==T)
		vals2keep = which(!is.na(cordCal_unix)==T)

		overlap_index = which(greg_unix %in% cordCal_unix)
		tas_365_minus1 = interp1(x=cordCal_unix[vals2keep], y=tas_all[vals2keep], xi = greg_unix[1:(length(greg_unix)-1)], method = c("linear"))

		# we miss very last dec 31 as x only goes to dec 30 and inter1 cannot interp beyond range of x (extrapolate). thats why we do this -1: xi = greg_unix[1:(length(greg_unix)-1)]
		# to fix we simply repeat last value
		tas_365 = c(tas_365_minus1, tas_365_minus1[length(tas_365_minus1)])

		#datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
		days_cordex=greg_cal_rcp
		years_rcp=as.numeric(substr(greg_cal_rcp,1,4))

		# do qmap
		rcp26 = doQmap(tas_365,pars)

		#annual aggregated timeseries
		rcp26_year=aggregate(rcp26, list(years_rcp), 'mean')

		# trim to common length
		si = which(greg_cal_rcp==startDateRcp26)
		ei = which(greg_cal_rcp==endDateRcp26)

		# assign to list
		rcp26_qmap_list[[rcp26_file]] <- list(rcp26[si:ei])
		rcp26_nqmap[[rcp26_file]] <- list(tas_365[si:ei])
		rcp26_dates[[rcp26_file]]<- list(greg_cal_rcp[si:ei])
		
		# do seasonal qmap
		rcp26_summer = doQmap(tas_365[summer_index],pars_summer)
		rcp26_autumn = doQmap(tas_365[autumn_index],pars_autumn)
		rcp26_winter= doQmap(tas_365[winter_index],pars_winter)
		rcp26_spring = doQmap(tas_365[spring_index],pars_spring)

		rcp26_season<-1:length(rcp26)
		rcp26_season[summer_index]<-rcp26_summer
		rcp26_season[autumn_index]<-rcp26_autumn
		rcp26_season[winter_index]<-rcp26_winter
		rcp26_season[spring_index]<-rcp26_spring


		rcp26_qmap_season_list[[rcp26_file]] <- list(rcp26_season[si:ei])
		}

		mydata=rcp26_qmap_list
		mydates=rcp26_dates
		modelNames = sapply(strsplit(names(mydata),'/'), "[", 9)
		df=as.data.frame(mydata)
		names(df)<-modelNames
		df$Date<-as.Date(unlist(mydates[[1]] ))
		save(df,file = paste0(outdir,"rcp26_",var))


	# Cordex rcp85 2006-2100=======================================================================
	rcp85_qmap_list <- list()
	rcp85_nqmap <- list()
	rcp85_dates <- list()
	rcp85_qmap_season_list <- list()
	for (rcp85_file in rcp85_files){

		modfile=rcp85_file
		# define variable and gridbox
		nc=nc_open(modfile)
		if(!'tas' %in% names(nc$var) ){print("not enough vars, skipping this file");print(length(nc$var));next}

		tas =ncvar_get(nc, 'tas') 		# !!HARDCODE!!
		lon =ncvar_get(nc, 'lon')
		lat =ncvar_get(nc, 'lat')
		myx =which.min(abs(lon - mylon))
		myy = which.min(abs(lat - mylat))

		# extract timeseries corresponding to qmap (all)
		tas_all = tas[myx,myy,]
		if(all(is.na(tas_all)==TRUE)){print("no data found, skipping this file");next}
		if(anyNA(tas_all)==TRUE){print("NAs found, skipping this file");next}

		# define time
		time = ncvar_get( nc,'time')
		time2=time
		z <- time2*24*60*60 #convert days to seconds
		#origin=substr(nc$dim$time$units,13,20)
		origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
		origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section

		# calender dependent section
		cal = nc$dim$time$calendar
		print(cal)

	# standard cal is assumed to equal proleptic_gregorian
		if(cal=='standard'){cal <- "proleptic_gregorian"}
		if(cal=='gregorian'){cal <- "proleptic_gregorian"}

		if(cal=="360_day"){
			cal <- "360_day"
			seconds.per.day <- 86400
			ts.dat.days <- time
			origin.pcict <- as.PCICt(origin, cal)
			ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
			datesPl <- ts.dat.pcict
			}
		if(cal=="365_day"){
			cal <- "365_day"
			seconds.per.day <- 86400
			ts.dat.days <- time
			origin.pcict <- as.PCICt(origin, cal)
			ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
			datesPl <- ts.dat.pcict
			}

		if(cal=="proleptic_gregorian"){
			datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
			}

		# convert 360 to 365/366 calender with interpolation
		#greg_cal = seq(as.Date('1970-01-01'), as.Date('2005-12-31'), "days") # hardcoded
		if(cal=="360_day"){
			greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(ENDDATE_CLIM), "days")
			}

		if(cal=="365_day"){
			greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") 
			}

		if(cal=="proleptic_gregorian"){
			greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(ENDDATE_CLIM), "days") # no +1 need
			}


		# get index for interp1d
		greg_unix = as.numeric(as.POSIXct(greg_cal_rcp, format="%Y-%m-%d"))
		cordCal_unix = as.numeric(as.POSIXct( as.Date(substr(datesPl,1,10)), format="%Y-%m-%d"))

		# NAs generated in cordCal_unix as dates are made in the pcict calender conversion that dont exists ie Feb 30
		# we remove these values completely representing 0.48% of values (63/12960)
		vals2rm = which(is.na(cordCal_unix)==T)
		vals2keep = which(!is.na(cordCal_unix)==T)


		overlap_index = which(greg_unix %in% cordCal_unix)
		tas_365_minus1 = interp1(x=cordCal_unix[vals2keep], y=tas_all[vals2keep], xi = greg_unix[1:(length(greg_unix)-1)], method = c("linear"))

		# we miss very last dec 31 as x only goes to dec 30 and inter1 cannot interp beyond range of x (extrapolate). thats why we do this -1: xi = greg_unix[1:(length(greg_unix)-1)]
		# to fix we simply repeat last value
		tas_365 = c(tas_365_minus1, tas_365_minus1[length(tas_365_minus1)])

		#datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
		days_cordex=greg_cal_rcp
		years_rcp=as.numeric(substr(greg_cal_rcp,1,4))


		# do qmap
		rcp85 = doQmap(tas_365,pars)

		#annual aggregated timeseries
		rcp85_year=aggregate(rcp85, list(years_rcp), 'mean')

			# trim to common length
		si = which(greg_cal_rcp==startDateRcp85)
		ei = which(greg_cal_rcp==endDateRcp85)

		# assign to list
		rcp85_qmap_list[[rcp85_file]] <- list(rcp85[si:ei])
		rcp85_nqmap[[rcp85_file]] <- list(tas_365[si:ei])
		rcp85_dates[[rcp85_file]]<- list(greg_cal_rcp[si:ei])

		# do seasonal qmap
		rcp85_summer = doQmap(tas_365[summer_index],pars_summer)
		rcp85_autumn = doQmap(tas_365[autumn_index],pars_autumn)
		rcp85_winter= doQmap(tas_365[winter_index],pars_winter)
		rcp85_spring = doQmap(tas_365[spring_index],pars_spring)

		rcp85_season<-1:length(rcp85)
		rcp85_season[summer_index]<-rcp85_summer
		rcp85_season[autumn_index]<-rcp85_autumn
		rcp85_season[winter_index]<-rcp85_winter
		rcp85_season[spring_index]<-rcp85_spring

		rcp85_qmap_season_list[[rcp85_file]] <- list(rcp85_season[si:ei])

		}

		mydata=rcp85_qmap_list
		mydates=rcp85_dates
		modelNames = sapply(strsplit(names(mydata),'/'), "[", 9)
		df=as.data.frame(mydata)
		names(df)<-modelNames
		df$Date<-as.Date(unlist(mydates[[1]] ))
		save(df,file = paste0(outdir,"rcp85_",var))


#===============================================================================
# Plots
#===============================================================================
mycol = wes_palette("Zissou1", 3, type = "continuous")

#===============================================================================
# prepare data yearly
#===============================================================================

# prepare data annual means and stats
df=as.data.frame(hist_qmap_list)
years_hist=( substr( as.Date(unlist(hist_dates[[1]] )), 1, 4)) 
hist_year = aggregate(df, list(years_hist), mean)
hist_year$SD=apply(hist_year[,2:dim(hist_year)[2]],1, sd, na.rm = TRUE)
hist_year$MEAN=apply(hist_year[,2:(dim(hist_year)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp26_qmap_list)
years_rcp26=( substr( as.Date(unlist(rcp26_dates[[1]] )), 1, 4)) 
rcp26_year = aggregate(df, list(years_rcp26), mean)
rcp26_year$SD=apply(rcp26_year[,2:dim(rcp26_year)[2]],1, sd, na.rm = TRUE)
rcp26_year$MEAN=apply(rcp26_year[,2:(dim(rcp26_year)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp85_qmap_list)
years_rcp85=( substr( as.Date(unlist(rcp85_dates[[1]] )), 1, 4)) 
rcp85_year = aggregate(df, list(years_rcp85), mean)
rcp85_year$SD=apply(rcp85_year[,2:dim(rcp85_year)[2]],1, sd, na.rm = TRUE)
rcp85_year$MEAN=apply(rcp85_year[,2:(dim(rcp85_year)[2]-1)],1, mean, na.rm = TRUE)

#===envelope plots========================================================================
if (var == 'tas'){

	pdf(paste(outdir,var,"TS2.pdf"))
	# caption: Mean annual near surface air temperature at the Weissfluhjoch (2540m asl)  showing corrected historical, RCP2.6 and RCP8.5 timseries. Observations and uncorrected historical data are also shown for comparison. The coloured envelops indicate +/- 1 SD of the model spread and multi-modal mean is given by the bold line.
	lwd=3
	plot(rcp26_year$Group.1, rcp26_year$MEAN, xlim=c(1970,2100),ylim=c(270,280), type='l', col=mycol[1], lwd=lwd, ylab="Air temperacture (K)", xlab=" ")


	y = c(rcp26_year$MEAN- rcp26_year$SD,rev(rcp26_year$MEAN+ rcp26_year$SD))
	x=c((rcp26_year$Group.1),rev(rcp26_year$Group.1))
	polygon (x,y, col=rgb(0, 0, 1,0.5),border='NA')



	lines(rcp85_year$Group.1, rcp85_year$MEAN, ylim=c(270,280), type='l', col=mycol[3],lwd=lwd)
	y = c(rcp85_year$MEAN- rcp85_year$SD,rev(rcp85_year$MEAN+ rcp85_year$SD))
	x=c((rcp85_year$Group.1),rev(rcp85_year$Group.1))
	polygon (x,y, col=rgb(1, 0, 0,0.5),border='NA')


	lines(hist_year$Group.1, hist_year$MEAN, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd)
	y = c(hist_year$MEAN- hist_year$SD,rev(hist_year$MEAN+ hist_year$SD))
	x=c((hist_year$Group.1),rev(hist_year$Group.1))
	polygon (x,y, col=rgb(0, 1, 0,0.5),border='NA')

	lines(tas_365_year, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd,lty=2)
	lines(obs_year,lwd=lwd)
	legend("bottomright", c("Hist", "RCP2.6", "RCP8.5", "Hist_uncorrected", "OBS"), col=c(rgb(0, 1, 0),mycol[1],mycol[3],mycol[2], "black"),lty=c(1,1,1,2,1),lwd=3)
	dev.off()

}
#===========================================================================================







# Evaluation historical (1980-01-01 - 2005-12-31)
# prepare data annual means and stats
#
#hist_noqmap = hist_cp
#hist_qmap =hist[start_cp_cord:end_cp_cord]


#hist_month = aggregate(hist_noqmap, list(dates_month), mean)[,2]
#hist_qmap_month = aggregate(hist_qmap, list(dates_month), mean)[,2]


#===============================================================================
# prepare data monthly
#===============================================================================

# prepare data monthly  means and stats historical for evaluation


dates_cp = greg_cal_hist[start_cp_cord:end_cp_cord]
dates_month = substring(dates_cp,1,7)

obs = obs_daily[1:length(hist_cp)]
obs_month = aggregate(obs, list(dates_month), mean)[,2]


# nned to redo indexs
si = which(as.Date(unlist(hist_dates[[1]]))==dates_cp[1])
ei = which(as.Date(unlist(hist_dates[[1]]))==dates_cp[length(dates_cp)])
df=as.data.frame(hist_nqmap)
hist_month = aggregate(df[si:ei,], list(dates_month), mean,na.rm = TRUE)
hist_month$SD=apply(hist_month[,2:dim(hist_month)[2]],1, sd, na.rm = TRUE)
hist_month$MEAN=apply(hist_month[,2:(dim(hist_month)[2]-1)],1, mean, na.rm = TRUE)



df=as.data.frame(hist_qmap_list)
hist_qmap_month = aggregate(df[si:ei,], list(dates_month), mean,na.rm = TRUE)
hist_qmap_month$SD=apply(hist_qmap_month[,2:dim(hist_qmap_month)[2]],1, sd, na.rm = TRUE)
hist_qmap_month$MEAN=apply(hist_qmap_month[,2:(dim(hist_qmap_month)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(hist_qmap_season_list)
hist_qmap_season_month = aggregate(df[si:ei,], list(dates_month), mean,na.rm = TRUE)
#hist_qmap_season_month$SD=apply(hist_qmap_season_month[,2:dim(hist_qmap_season_month)[2]],1, sd, na.rm = TRUE)
#hist_qmap_season_month$MEAN=apply(hist_qmap_season_month[,2:(dim(hist_qmap_season_month)[2]-1)],1, mean, na.rm = TRUE)

#df=as.data.frame(rcp26_qmap_list)
#rcp26_year = aggregate(df, list(years_rcp), mean)
#rcp26_year$SD=apply(rcp26_year[,2:dim(rcp26_year)[2]],1, sd, na.rm = TRUE)
#rcp26_year$MEAN=apply(rcp26_year[,2:(dim(rcp26_year)[2]-1)],1, mean, na.rm = TRUE)

#df=as.data.frame(rcp85_qmap_list)
#rcp85_year = aggregate(df, list(years_rcp), mean)
#rcp85_year$SD=apply(rcp85_year[,2:dim(rcp85_year)[2]],1, sd, na.rm = TRUE)
#rcp85_year$MEAN=apply(rcp85_year[,2:(dim(rcp85_year)[2]-1)],1, mean, na.rm = TRUE)



#===============================================================================
# summary stats
#===============================================================================
cor(obs_month, hist_month$MEAN)
cor(obs_month, hist_qmap_month$MEAN)

rmse(obs_month, hist_month$MEAN)
rmse(obs_month, hist_qmap_month$MEAN)


mean(obs_month)
mean(hist_month$MEAN)
mean(hist_qmap_month$MEAN)

# plot(obs_month, hist_month$MEAN)
# plot(obs_month, hist_qmap_month$MEAN)

#===============================================================================
# ecdf
#===============================================================================
lwd=3
# plot(obs_month,type='l', ylim = c(255,290), lwd=lwd)
# lines( hist_month$MEAN,type='l', col='red', lwd=lwd)
# lines(hist_qmap_month$MEAN,type='l', col='blue', lwd=lwd)
#legend("bottomright", c("OBS", "CLIM", "CLIM_QM"), col=c("black","red", "blue"), lwd=lwd)
#pdf(paste(outdir,var,"CDF2.pdf"))
 plot(ecdf(obs_month), ylab="cdf",  col=mycol[1], xlab="TA (K)", main=var, lwd=lwd,lty=1)
 # lines(ecdf(hist_month$MEAN),  col=mycol[2], lwd=lwd)
 # lines(ecdf(hist_qmap_month$MEAN),  col=mycol[3], lwd=1)

  lines(ecdf(hist_month[,2]),  col=mycol[2], lwd=lwd)
 lines(ecdf(hist_qmap_month[,2]),  col=mycol[3], lwd=1)
 lines(ecdf(hist_qmap_season_month[,2]),  col='green', lwd=1)
legend("topright", c("OBS", "CLIM", "CLIM_QM"), col=c(mycol[1],mycol[2], mycol[3]),lty=c(1,1,1) ,lwd=lwd)
#dev.off()
# plot(ecdf(obs), ylab="cdf", xlab="TA (K)", main=" ", lwd=lwd)
# lines(ecdf(hist_noqmap), col='red', lwd=lwd)
# lines(ecdf(hist_qmap), col='blue', lwd=lwd)
#legend("bottomright", c("OBS", "CLIM", "CLIM_QM"), col=c("black","red", "blue"), lwd=lwd)


#===============================================================================
# doy
#===============================================================================

# prepare daily data
hist=as.data.frame(hist_nqmap)
hist$SD=apply(hist[,2:dim(hist)[2]],1, sd, na.rm = TRUE)
hist$MEAN=apply(hist[,2:(dim(hist)[2]-1)],1, mean, na.rm = TRUE)

hist_qmap=as.data.frame(hist_qmap_list)
hist_qmap$SD=apply(hist_qmap[,2:dim(hist_qmap)[2]],1, sd, na.rm = TRUE)
hist_qmap$MEAN=apply(hist_qmap[,2:(dim(hist_qmap)[2]-1)],1, mean, na.rm = TRUE)

hist_qmap_season=as.data.frame(hist_qmap_season_list)
# make doy vector
doy=substr(unique(dates_cp), 6,10)


# aggeragte to DOY and cut to cp common period
doy_obs = aggregate(obs, list(doy), mean)
doy_clim = aggregate(hist$MEAN[si:ei], list(doy), mean)
doy_qm = aggregate(hist_qmap$MEAN[si:ei], list(doy), mean)

doy_clim = aggregate(hist[si:ei,2], list(doy), mean)
doy_qm = aggregate(hist_qmap[si:ei,2], list(doy), mean)
doy_qm_season = aggregate(hist_qmap_season[si:ei,2], list(doy), mean)
#
#pdf(paste(outdir,var,"DOY2.pdf"))
plot(doy_obs$x, type='l',col=mycol[1],lwd=lwd,lty=2, xlab='DOY', ylab='Air temperature [k]', main=var)
lines(doy_clim$x, type='l', col=mycol[2], lwd=lwd)
lines(doy_qm$x, type='l', col=mycol[3],lwd=lwd)
lines(doy_qm_season$x, type='l', col='green',lwd=lwd)
legend("topright", c("OBS", "CLIM", "CLIM_QM"), col=c(mycol[1],mycol[2], mycol[3]),lty=c(2,1,1) ,lwd=lwd)
#dev.off()

#cor plots
# par(mfrow=c(1,2))
# plot(doy_obs$x,doy_clim$x,col=mycol[1])
# rmse(doy_obs$x,doy_clim$x)
# abline(0,1)
# plot(doy_obs$x,doy_qm$x,col=mycol[1])
# rmse(doy_obs$x,doy_qm$x)
# abline(0,1)




# we evaluate with rcp data in period 2006-2017
}
dev.off()
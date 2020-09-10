require(qmap)
require(ncdf4)

require(caTools)
require(PCICt) # converts to 365 day clim model calender
require(pracma)
require(hydroGOF)
require(wesanderson)

# create 1h data by
# import pandas as pd
# path_inp = "/home/joel/sim/qmap/wfj_long.csv"
# # lesson of psum mess is dont use psum!!
# df_obs= pd.read_csv(path_inp, index_col=0, parse_dates=True)
# df.resample('1H').interpolate()
# df_1h = df_obs.resample('1H').interpolate()
# df_1hto_csv(path_or_buf='/home/joel/sim/qmap/wfj_long_1H.csv' ,na_rep=-999,float_format='%.8f', header=True, sep=',')
		
# TAS script
# handle different calender
# Setup ========================================================================

# POI (WFJ) - to extract cordex gridbox
mylon = 9.80490
mylat = 46.83066

# set end date of climate timeseires to cut to as different model provide differnet length timeseries
ENDDATE_CLIM="2099-12-31"

# files
indir="/home/joel/sim/qmap/topoclim"
outdir=paste0(indir,"/fsm2/")
files = list.files(path=indir, pattern="HOURLY.txt$", recursive=T, full.name=T)
dir.create(outdir)
hist_files = files[ grep('historical',files)]
rcp26_files = files[ grep('rcp26',files)]
rcp85_files = files[ grep('rcp85',files)]

# define periods
startDateHist = "1979-01-01 00:00:00"
endDateHist = "2005-12-31 23:00:00"
startDateRcp26 = "2006-01-02 00:00:00"
endDateRcp26 = "2099-12-31 23:00:00"
startDateRcp85 = "2006-01-02 00:00:00"
endDateRcp85 = "2099-12-31 23:00:00"

# Cordex data
#hist_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/historical/r1i1p1/REMO2015/v1/day/tas/tas_EUR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"
#rcp26_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/rcp26/r1i1p1/REMO2015/v1/day/tas/tas_EUR-22_MOHC-HadGEM2-ES_rcp26_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"
#rcp85_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/rcp85/r1i1p1/REMO2015/v1/day/tas/tas_EUR-22_MOHC-HadGEM2-ES_rcp85_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"



# Obs Era5 data downscaled by toposcale, can be resampled
obsfile = "/home/joel/sim/qmap/wfj_long_1H.csv"#
obs=read.csv(obsfile)
myvars=c('TA', 'RH', 'VW', 'DW', 'P', 'ISWR', 'ILWR','PINT')#,  'pr', 'uas', 'vas')







	# Cordex historical "1970-01-01 12" UTC" to 2005-12-30 12"======================

# lists of parameters to be populate
pars_ta<-list()
pars_rh<-list()
pars_vw<-list()
pars_dw<-list()
pars_p<-list()
pars_iswr<-list()
pars_ilwr<-list()
pars_pint<-list()

pars_summer_ta<-list()
pars_summer_rh<-list()
pars_summer_vw<-list()
pars_summer_dw<-list()
pars_summer_p<-list()
pars_summer_iswr<-list()
pars_summer_ilwr<-list()
pars_summer_pint<-list()

pars_autumn_ta<-list()
pars_autumn_rh<-list()
pars_autumn_vw<-list()
pars_autumn_dw<-list()
pars_autumn_P<-list()
pars_autumn_iswr<-list()
pars_autumn_ilwr<-list()
pars_autumn_pint<-list()

pars_winter_ta<-list()
pars_winter_rh<-list()
pars_winter_vw<-list()
pars_winter_dw<-list()
pars_winter_P<-list()
pars_winter_iswr<-list()
pars_winter_ilwr<-list()
pars_winter_pint<-list()

pars_spring_ta<-list()
pars_spring_rh<-list()
pars_spring_vw<-list()
pars_spring_dw<-list()
pars_spring_P<-list()
pars_spring_iswr<-list()
pars_spring_ilwr<-list()
pars_spring_pint<-list()


for (hist_file in hist_files){
	stop = FALSE
	print(hist_file)
	modelNameBase = unlist(strsplit(hist_file,'/'))[length(unlist(strsplit(hist_file,'/')))]
	GCM = unlist(strsplit(modelNameBase,'_historical'))[1]
	RCM = unlist(strsplit(modelNameBase,'_'))[4]
	modelChain=paste0(GCM,'_',RCM)

	nc=read.csv(hist_file, sep=',', header=T, na.strings = "-999")

	hist_qmap_list <- list()
	hist_nqmap <- list()
	hist_dates<-list()
	hist_qmap_season_list <- list()

	rcp26_qmap_list <- list()
	rcp26_nqmap <- list()
	rcp26_dates <-list()
	rcp26_qmap_season_list <- list()

	rcp85_qmap_list <- list()
	rcp85_nqmap <- list()
	rcp85_dates <-list()
	rcp85_qmap_season_list <- list()

	for (var in myvars){
	print(var)
	# Obs 1980-01-01 - 2017-12-31 ==========================================================================
	#read obs variable TA


		if(var=="TA"){obsindex <-11; convFact <-1; modindex <- 2} # K
			if(var=="RH"){obsindex <-8; convFact <-100; modindex <- 3} # kgm-2s-1
				if(var=="VW"){obsindex <-12; convFact <-1; modindex <- 4}# Pa
					if(var=="DW"){obsindex <-2; convFact <-1; modindex <- 5} # % 0-100
						if(var=="P"){obsindex <-5; convFact <-1; modindex <- 6}	# Wm-2
							if(var=="ISWR"){obsindex <-4; convFact <-1; modindex <- 7}# Wm-2
								if(var=="ILWR"){obsindex <-3; convFact <-1; modindex <- 8}# ms-1
									if(var=="PINT"){obsindex <-6; convFact <-1; modindex <- 9}# ms-1

		obs_var=obs[,obsindex]*convFact	# !!HARDCODE!!

		# aggregate obs to daily resolution
		obs_datetime= strptime(obs$datetime,format="%Y-%m-%d %H:%M:%S")


		tas =nc[,modindex]

		if (sum(is.na(tas))==length(tas)){
			print ("No data found, skipping file");next
		}
		


		tas_all = tas
		if(all(is.na(tas_all)==TRUE)){print("no data found, skipping this file");stop=TRUE;next}

		#if length(which(is.na(tas)))  > 50
		if(anyNA(tas_all)==TRUE){print("NAs found, skipping this file");stop=TRUE;next}
		# define time

		datesPl= nc[,1]
		greg_cal_hist = strptime( datesPl,format="%Y-%m-%d %H:%M:%S")#seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") 

		days_cordex=greg_cal_hist
		years_hist=as.numeric(substr(greg_cal_hist,1,4))
		start_cp_cord = which(strptime(days_cordex, format="%Y-%m-%d %H:%M:%S")==(obs_datetime[1]))
		end_cp_cord=which( strptime(days_cordex, format="%Y-%m-%d %H:%M:%S")==(obs_datetime[length(obs_datetime)])   )

		if(length(end_cp_cord)==0){ # if true this mean that obs extend beyond end_cp_cord of the cordex data - this is the case for cordex historical
			end_cp_cord= length(days_cordex) # then simply take last cordex date
			}
		if(length(start_cp_cord)==0){ # if true this mean that obs starts before cordex data - this is the case for cordex rcp data
			start_cp_cord = 1 # then simply take first cordex date
			}
		# extract timeseries corresponding to obs (subset in time/space)

		hist_cp=tas[start_cp_cord:end_cp_cord]

		# define common period
		start_cp= greg_cal_hist[start_cp_cord]
		end_cp = greg_cal_hist[end_cp_cord]

		# cut obs to common period cp
		start_obs = which(obs_datetime==start_cp)
		end_obs=which(obs_datetime==end_cp)
		obs_daily_cp = obs_var[start_obs: end_obs]

		# get gmap pars
		pars = fitQmap(obs_daily_cp,hist_cp, method = "QUANT")

		# do qmap
		hist = doQmap(tas,pars)

		#annual aggregated timeseries
		hist_year=aggregate(hist, list(years_hist), 'mean')
		tas_year=aggregate(tas, list(years_hist), 'mean')
		obs_year=aggregate(obs_var, list(substring(obs_datetime,1,4)), 'mean')

		# trim to common length
		si = which(greg_cal_hist==startDateHist)
		ei = which(greg_cal_hist==endDateHist)

		# add to list
		hist_qmap_list[[var]] <- list(hist[si:ei])
		hist_nqmap[[var]] <- list(tas[si:ei])
		hist_dates[[var]]<- list(greg_cal_hist[si:ei] )
		



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

		
		hist_summer = doQmap(tas[summer_index],pars_summer)
		hist_autumn = doQmap(tas[autumn_index],pars_autumn)
		hist_winter= doQmap(tas[winter_index],pars_winter)
		hist_spring = doQmap(tas[spring_index],pars_spring)

		hist_season<-1:length(hist)
		hist_season[summer_index]<-hist_summer
		hist_season[autumn_index]<-hist_autumn
		hist_season[winter_index]<-hist_winter
		hist_season[spring_index]<-hist_spring


		hist_qmap_season_list[[var]] <- list(hist_season[si:ei])
		

}




		if (stop){next}
		modelName = unlist(strsplit(hist_file,'/'))[ 7]
		df=as.data.frame(hist_qmap_season_list)
		names(df)<-names(hist_qmap_season_list)
		datetime<-greg_cal_hist[si:ei] 
		df2=data.frame(datetime, round(df[1:7],1), round(df[8],4))
		write.csv(df2,file = paste0(outdir,modelName, "_QMAP.txt"), row.names=FALSE)



}

	# Cordex rcp26 =======================================================================


	for (rcp26_file in rcp26_files){
		
	stop = FALSE
	print(rcp26_file)
	nc=read.csv(rcp26_file, sep=',', header=T, na.strings = "-999")

	rcp26_qmap_list <- list()
	rcp26_nqmap <- list()
	rcp26_dates <-list()
	rcp26_qmap_season_list <- list()

	for (var in myvars){
		print(var)
		if(var=="TA"){ modindex <- 2} # K
			if(var=="RH"){modindex <- 3} # kgm-2s-1
				if(var=="VW"){modindex <- 4}# Pa
					if(var=="DW"){modindex <- 5} # % 0-100
						if(var=="P"){modindex <- 6}	# Wm-2
							if(var=="ISWR"){modindex <- 7}# Wm-2
								if(var=="ILWR"){modindex <- 8}# ms-1
									if(var=="PINT"){ modindex <- 9}# ms-1

		tas =nc[,modindex]

		#print(head(tas))
		if (sum(is.na(tas))==length(tas)){print ("No data found, skipping file");stop=TRUE;next }

		tas_all = tas
		if(all(is.na(tas_all)==TRUE)){print("no data found, skipping this file");stop=TRUE;next}
		if(anyNA(tas_all)==TRUE){print("NAs found, skipping this file");stop=TRUE;next}
		# define time

		datesPl= nc[,1]
		greg_cal_rcp = strptime( datesPl,format="%Y-%m-%d %H:%M:%S")#seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") 

		days_cordex=greg_cal_rcp
		years_rcp=as.numeric(substr(greg_cal_rcp,1,4))

		# do qmap
		rcp26 = doQmap(tas,pars)

		#annual aggregated timeseries
		rcp26_year=aggregate(rcp26, list(years_rcp), 'mean')

		# trim to common length
		si = which(greg_cal_rcp==startDateRcp26)
		ei = which(greg_cal_rcp==endDateRcp26)

		# assign to list
		rcp26_qmap_list[[var]] <- list(rcp26[si:ei])
		rcp26_nqmap[[var]] <- list(tas[si:ei])
		rcp26_dates[[var]]<- list(greg_cal_rcp[si:ei])
		
		# do seasonal qmap
		rcp26_summer = doQmap(tas[summer_index],pars_summer)
		rcp26_autumn = doQmap(tas[autumn_index],pars_autumn)
		rcp26_winter= doQmap(tas[winter_index],pars_winter)
		rcp26_spring = doQmap(tas[spring_index],pars_spring)

		rcp26_season<-1:length(rcp26)
		rcp26_season[summer_index]<-rcp26_summer
		rcp26_season[autumn_index]<-rcp26_autumn
		rcp26_season[winter_index]<-rcp26_winter
		rcp26_season[spring_index]<-rcp26_spring


		rcp26_qmap_season_list[[var]] <- list(rcp26_season[si:ei])
		}
		


		if (stop){next}
		modelName = unlist(strsplit(rcp26_file,'/'))[ 7]
		df=as.data.frame(rcp26_qmap_season_list)
		names(df)<-names(rcp26_qmap_season_list)
		datetime<-greg_cal_rcp[si:ei]
		df2=data.frame(datetime, round(df[1:7],1), round(df[8],4))
		write.csv(df2,file = paste0(outdir,modelName, "_QMAP.txt"), row.names=FALSE)
	}


		# Cordex rcp26 =======================================================================


	for (rcp85_file in rcp85_files){
		
	stop = FALSE
	print(rcp85_file)
	nc=read.csv(rcp85_file, sep=',', header=T, na.strings = "-999")

	rcp85_qmap_list <- list()
	rcp85_nqmap <- list()
	rcp85_dates <-list()
	rcp85_qmap_season_list <- list()

	for (var in myvars){
		print(var)
		if(var=="TA"){ modindex <- 2} # K
			if(var=="RH"){modindex <- 3} # kgm-2s-1
				if(var=="VW"){modindex <- 4}# Pa
					if(var=="DW"){modindex <- 5} # % 0-100
						if(var=="P"){modindex <- 6}	# Wm-2
							if(var=="ISWR"){modindex <- 7}# Wm-2
								if(var=="ILWR"){modindex <- 8}# ms-1
									if(var=="PINT"){ modindex <- 9}# ms-1

		tas =nc[,modindex]

		#print(head(tas))
		if (sum(is.na(tas))==length(tas)){print ("No data found, skipping file");stop=TRUE;next }

		tas_all = tas
		if(all(is.na(tas_all)==TRUE)){print("no data found, skipping this file");stop=TRUE;next}
		if(anyNA(tas_all)==TRUE){print("NAs found, skipping this file");stop=TRUE;next}
		# define time

		datesPl= nc[,1]
		greg_cal_rcp = strptime( datesPl,format="%Y-%m-%d %H:%M:%S")#seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10)), "days") 

		days_cordex=greg_cal_rcp
		years_rcp=as.numeric(substr(greg_cal_rcp,1,4))

		# do qmap
		rcp85 = doQmap(tas,pars)

		#annual aggregated timeseries
		rcp85_year=aggregate(rcp85, list(years_rcp), 'mean')

		# trim to common length
		si = which(greg_cal_rcp==startDateRcp26)
		ei = which(greg_cal_rcp==endDateRcp26)

		# assign to list
		rcp85_qmap_list[[var]] <- list(rcp85[si:ei])
		rcp85_nqmap[[var]] <- list(tas[si:ei])
		rcp85_dates[[var]]<- list(greg_cal_rcp[si:ei])
		
		# do seasonal qmap
		rcp85_summer = doQmap(tas[summer_index],pars_summer)
		rcp85_autumn = doQmap(tas[autumn_index],pars_autumn)
		rcp85_winter= doQmap(tas[winter_index],pars_winter)
		rcp85_spring = doQmap(tas[spring_index],pars_spring)

		rcp85_season<-1:length(rcp85)
		rcp85_season[summer_index]<-rcp85_summer
		rcp85_season[autumn_index]<-rcp85_autumn
		rcp85_season[winter_index]<-rcp85_winter
		rcp85_season[spring_index]<-rcp85_spring


		rcp85_qmap_season_list[[var]] <- list(rcp85_season[si:ei])
		}



		if (stop){next}
		modelName = unlist(strsplit(rcp85_file,'/'))[ 7]
		df=as.data.frame(rcp85_qmap_season_list)
		names(df)<-names(rcp85_qmap_season_list)
		datetime<-greg_cal_rcp[si:ei]
		df2=data.frame(datetime, round(df[1:7],1), round(df[8],4))
		write.csv(df2,file = paste0(outdir,modelName, "_QMAP.txt"), row.names=FALSE)
	}






	
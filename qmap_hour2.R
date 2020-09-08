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
outdir=paste0(indir,"/aqmap_results/")
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

pdf(paste(outdir,"evalplot_singleModel_season.pdf"),height=20, width=10)
par(mfrow=c(8,2))

# Obs Era5 data downscaled by toposcale, can be resampled
obsfile = "/home/joel/sim/qmap/wfj_long_1H.csv"#
obs=read.csv(obsfile)
myvars=c('TA', 'RH', 'VW', 'DW', 'P', 'ISWR', 'ILWR','PINT')#,  'pr', 'uas', 'vas')







	# Cordex historical "1970-01-01 12" UTC" to 2005-12-30 12"======================


for (hist_file in hist_files){
	stop = FALSE
	print(hist_file)
	nc=read.csv(hist_file, sep=',', header=T, na.strings = "-999")

	hist_qmap_list <- list()
	hist_nqmap <- list()
	hist_dates<-list()
	hist_qmap_season_list <- list()

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



#===============================================================================
# Plots
#===============================================================================
mycol = wes_palette("Zissou1", 3, type = "continuous")

#===============================================================================
# prepare data yearly
#===============================================================================

# prepare data annual means and stats
df=as.data.frame(hist_qmap_list)
years_hist=( substr(hist_dates, 1, 4)) 
hist_year = aggregate(df, list(years_hist), mean)
hist_year$SD=apply(hist_year[,2:dim(hist_year)[2]],1, sd, na.rm = TRUE)
hist_year$MEAN=apply(hist_year[,2:(dim(hist_year)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp26_qmap_list)

#greg_cal_rcp[si:ei]

years_rcp26=( substr( greg_cal_rcp[si:ei], 1, 4)) 
rcp26_year = aggregate(df, list(years_rcp26), mean)
rcp26_year$SD=apply(rcp26_year[,2:dim(rcp26_year)[2]],1, sd, na.rm = TRUE)
rcp26_year$MEAN=apply(rcp26_year[,2:(dim(rcp26_year)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp85_qmap_list)
years_rcp85=( substr( greg_cal_rcp[si:ei], 1, 4)) 
rcp85_year = aggregate(df, list(years_rcp85), mean)
rcp85_year$SD=apply(rcp85_year[,2:dim(rcp85_year)[2]],1, sd, na.rm = TRUE)
rcp85_year$MEAN=apply(rcp85_year[,2:(dim(rcp85_year)[2]-1)],1, mean, na.rm = TRUE)

#===envelope plots========================================================================
if (var == 'TA'){

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

	lines(tas_year, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd,lty=2)
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

obs = obs_var[1:length(hist_cp)]
obs_month = aggregate(obs, list(dates_month), mean)[,2]


# nned to redo indexs
si = which(hist_dates==dates_cp[1])
ei = which(hist_dates==dates_cp[length(dates_cp)])
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
legend("topright", c("OBS", "CLIM", "CLIM_QM", "CLIM_QMSEASON"), col=c(mycol[1],mycol[2], mycol[3], 'green'),lty=c(1,1,1,1) ,lwd=lwd)
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

# redefine to fit length why needed??
obs = obs_var[1:(length(hist_cp)-25)]
# aggeragte to DOY and cut to cp common period
doy_obs = aggregate(obs, list(doy), mean)
doy_clim = aggregate(hist$MEAN[si:(ei-25)], list(doy), mean)
doy_qm = aggregate(hist_qmap$MEAN[si:(ei-25)], list(doy), mean)

doy_clim = aggregate(hist[si:(ei-25),2], list(doy), mean)
doy_qm = aggregate(hist_qmap[si:(ei-25),2], list(doy), mean)
doy_qm_season = aggregate(hist_qmap_season[si:(ei-25),2], list(doy), mean)
#
#pdf(paste(outdir,var,"DOY2.pdf"))
plot(doy_obs$x, type='l',col=mycol[1],lwd=lwd,lty=2, xlab='DOY', ylab='Air temperature [k]', main=var)
lines(doy_clim$x, type='l', col=mycol[2], lwd=lwd)
lines(doy_qm$x, type='l', col=mycol[3],lwd=lwd)
lines(doy_qm_season$x, type='l', col='green',lwd=lwd)
#legend("topright", c("OBS", "CLIM", "CLIM_QM"), col=c(mycol[1],mycol[2], mycol[3]),lty=c(2,1,1) ,lwd=lwd)
legend("topright", c("OBS", "CLIM", "CLIM_QM", "CLIM_QMSEASON"), col=c(mycol[1],mycol[2], mycol[3], 'green'),lty=c(2,1,1,1) ,lwd=lwd)
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

#=====================================================================
# generate sim files
#=====================================================================
		 			# 	"year": dates.year, 
		 			# 	"month": dates.month, 
						# "day": dates.day, 
						# "hour": dates.hour,
						# "ISWR":station.data_disagg.glob, 
						# "ILWR":station.data_disagg.ILWR, 
						# "Sf":snowTot/(60.*60.), # prate in mm/hr to kgm2/s
						# "Rf":rainTot/(60.*60.), # prate in mm/hr to kgm2/s
						# "TA":station.data_disagg.temp, 
						# "RH":station.data_disagg.hum,#*0.01, #meteoio 0-1
						# "VW":station.data_disagg.wind,
						# "P":station.data_disagg.P,

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


	
dir.create(paste0(outdir,  "/fsm_met/"))
hist_files = list.files(path=outdir, pattern="hist_", recursive=T, full.name=T)

load(paste0(outdir, "/hist_TA"))
TA=df
load(paste0(outdir, "/hist_RH"))
RH=df
load(paste0(outdir, "/hist_P"))
P=df
load(paste0(outdir, "/hist_VW"))
VW=df
load(paste0(outdir, "/hist_DW"))
DW=df
load(paste0(outdir, "/hist_ILWR"))
ILWR=df
load(paste0(outdir, "/hist_ISWR"))
ISWR=df

# for some reason ISWR has non zero min so rebase this
ISWR=ISWR-min(ISWR)
load(paste0(outdir, "/hist_PINT"))
PINT=df

histdate = TA$Date
year = substr(histdate, 1,4)  
month = substr(histdate, 6,7)  
day = substr(histdate, 9,10)  
hour  = substr(histdate, 12,13)  
models = names(TA)

for (i in 1:(dim(TA)[2]-1)){
TA_col <- TA[,i]
RH_col <- RH[,i]
P_col <- P[,i]
VW_col <- VW[,i]
DW_col <- DW[,i]
ILWR_col <- ILWR[,i]
ISWR_col <- ISWR[,i]
PINT_col <- PINT[,i]	

df=data.frame(histdate, TA_col, RH_col, P_col, VW_col, DW_col, ILWR_col, ISWR_col,PINT_col)
write.table(df, paste0(outdir, "/fsm_met/", models[i], "FSM.txt"), sep="\t",row.names=F, col.names=F)
}





load(paste0(outdir, "/rcp26_TA"))
TA=df
load(paste0(outdir, "/rcp26_RH"))
RH=df
load(paste0(outdir, "/rcp26_P"))
P=df
load(paste0(outdir, "/rcp26_VW"))
VW=df
load(paste0(outdir, "/rcp26_DW"))
DW=df
load(paste0(outdir, "/rcp26_ILWR"))
ILWR=df
load(paste0(outdir, "/rcp26_ISWR"))
ISWR=df
load(paste0(outdir, "/rcp26_PINT"))
PINT=df

rcpdate = TA$Date
year = substr(rcpdate, 1,4)  
month = substr(rcpdate, 6,7)  
day = substr(rcpdate, 9,10)  
hour  = substr(rcpdate, 12,13)  
models = names(TA)

for (i in 1:(dim(TA)[2]-1)){
TA_col <- TA[,i]
RH_col <- RH[,i]
P_col <- P[,i]
VW_col <- VW[,i]
DW_col <- DW[,i]
ILWR_col <- ILWR[,i]
ISWR_col <- ISWR[,i]
PINT_col <- PINT[,i]	

df=data.frame(rcpdate, TA_col, RH_col, P_col, VW_col, DW_col, ILWR_col, ISWR_col,PINT_col)
write.table(df, paste0(outdir, "/fsm_met/", models[i], "FSM.txt"), sep="\t",row.names=F, col.names=F)
}

load(paste0(outdir, "/rcp85_TA"))
TA=df
load(paste0(outdir, "/rcp85_RH"))
RH=df
load(paste0(outdir, "/rcp85_P"))
P=df
load(paste0(outdir, "/rcp85_VW"))
VW=df
load(paste0(outdir, "/rcp85_DW"))
DW=df
load(paste0(outdir, "/rcp85_ILWR"))
ILWR=df
load(paste0(outdir, "/rcp85_ISWR"))
ISWR=df
load(paste0(outdir, "/rcp85_PINT"))
PINT=df

rcpdate = TA$Date
year = substr(rcpdate, 1,4)  
month = substr(rcpdate, 6,7)  
day = substr(rcpdate, 9,10)  
hour  = substr(rcpdate, 12,13)  
models = names(TA)

for (i in 1:(dim(TA)[2]-1)){
TA_col <- TA[,i]
RH_col <- RH[,i]
P_col <- P[,i]
VW_col <- VW[,i]
DW_col <- DW[,i]
ILWR_col <- ILWR[,i]
ISWR_col <- ISWR[,i]
PINT_col <- PINT[,i]	

df=data.frame(year, month, day, hour, ISWR_col, ILWR_col, TA_col, RH_col, P_col, VW_col, PINT_col)
write.table(df, paste0(outdir, "/fsm_met/", models[i], "FSM.txt"), sep="\t",row.names=F, col.names=F)
}

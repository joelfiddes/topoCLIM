require(qmap)
require(ncdf4)
require(caTools)
require(pracma)
require(hydroGOF)
require(wesanderson)

# create 1h wfj era5 obs data by
# import pandas as pd
# path_inp = "/home/joel/sim/qmap/wfj_long.csv"
# # lesson of psum mess is dont use psum!!
# df_obs= pd.read_csv(path_inp, index_col=0, parse_dates=True)
# df.resample('1H').interpolate()
# df_1h = df_obs.resample('1H').interpolate()
# df_1hto_csv(path_or_buf='/home/joel/sim/qmap/wfj_long_1H.csv' ,na_rep=-999,float_format='%.8f', header=True, sep=',')
		
plot=FALSE
# Setup ========================================================================

# POI (WFJ) - to extract cordex gridbox
mylon = 9.80490
mylat = 46.83066


# directories
indir="/home/joel/sim/qmap/topoclim_test"
outdir=paste0(indir,"/aqmap_results/")
dir.create(outdir)

# dissaggregated calender corrected hourly input files (produced by topoCLIM.py)
files = list.files(path=indir, pattern="HOURLY.txt$", recursive=T, full.name=T)
hist_files = files[ grep('historical',files)]
rcp26_files = files[ grep('rcp26',files)]
rcp85_files = files[ grep('rcp85',files)]

# different models have different time spans so need to define common periods where all models have data

# common historical period of models
startDateHist = "1979-01-01 00:00:00"
endDateHist = "2005-12-31 23:00:00"

# common rcp periods of models - perhaps these can/should be 1 set, in practice is the same
startDateRcp26 = "2006-01-02 00:00:00"
endDateRcp26 = "2099-12-31 23:00:00"
startDateRcp85 = "2006-01-02 00:00:00"
endDateRcp85 = "2099-12-31 23:00:00"

# the period to do quantile mapping over (obs and sim must overlap!)
startDateQmap = "1990-01-01 00:00:00"
endDateQmap = "1999-12-31 23:00:00"


pdf(paste(outdir,"evalplot_singleModel_season.pdf"),height=20, width=10)
par(mfrow=c(8,2))

# Obs Era5 data downscaled by toposcale, can be resampled
obsfile = "/home/joel/sim/qmap/wfj_long_1H.csv"#
obs=read.csv(obsfile)
myvars=c('TA', 'RH', 'VW', 'DW', 'P', 'ISWR', 'ILWR','PINT')

for (var in myvars){
print(var)

	if(var=="TA"){obsindex <-11; convFact <-1; modindex <- 2} # K
		if(var=="RH"){obsindex <-8; convFact <-100; modindex <- 3} # 0-1 -> 0-100
			if(var=="VW"){obsindex <-12; convFact <-1; modindex <- 4}# ms-1
				if(var=="DW"){obsindex <-2; convFact <-1; modindex <- 5} # ms-1
					if(var=="P"){obsindex <-5; convFact <-1; modindex <- 6}	# Pa
						if(var=="ISWR"){obsindex <-4; convFact <-1; modindex <- 7}# Wm-2
							if(var=="ILWR"){obsindex <-3; convFact <-1; modindex <- 8}# Wm-2
								if(var=="PINT"){obsindex <-6; convFact <-1; modindex <- 9}# mm/hr
	# read and convert obs unit
	obs_var=obs[,obsindex]*convFact	

	# aggregate obs to daily resolution
	obs_datetime= strptime(obs$datetime,format="%Y-%m-%d %H:%M:%S")

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


		print(hist_file)
		modelNameBase = unlist(strsplit(hist_file,'/'))[length(unlist(strsplit(hist_file,'/')))]
		GCM = unlist(strsplit(modelNameBase,'_historical'))[1]
		RCM = unlist(strsplit(modelNameBase,'_'))[4]
		modelChain=paste0(GCM,'_',RCM)
		
		
		# get data
		nc=read.csv(hist_file, sep=',', header=T, na.strings = "-999")
		tas =nc[,modindex]

		if(all(is.na(tas)==TRUE)){print(paste0("no data found in variable, skipping this file", hist_file));next}
		
		# if all present and correct can now add to chain
		modelChain_vec = c(modelChain_vec,modelChain) # available hist data

		# define periods
		cordex_dates = strptime( nc[,1],format="%Y-%m-%d %H:%M:%S")
		start_cp_cord = which(strptime(cordex_dates, format="%Y-%m-%d %H:%M:%S")==startDateHist)
		end_cp_cord=which( strptime(cordex_dates, format="%Y-%m-%d %H:%M:%S")==endDateHist  )
		start_qmap_cord = which(strptime(cordex_dates, format="%Y-%m-%d %H:%M:%S")==startDateQmap)
		end_qmap_cord=which( strptime(cordex_dates, format="%Y-%m-%d %H:%M:%S")==endDateQmap  )

		
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

		hist_season<-1:length(hist)
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
	save(df,file = paste0(outdir,"hist_",var))

	# Cordex rcp26 =======================================================================
	rcp26_qmap_list <- list()
	rcp26_nqmap <- list()
	rcp26_qmap_season_list <- list()

	for (rcp26_file in rcp26_files){

		print(rcp26_file)
		modelNameBase = unlist(strsplit(rcp26_file,'/'))[length(unlist(strsplit(rcp26_file,'/')))]
		GCM = unlist(strsplit(modelNameBase,'_rcp26_'))[1]
		RCM = unlist(strsplit(modelNameBase,'_'))[4]
		modelChain=paste0(GCM,'_',RCM)

		# skip to next model if does not have hist data
		if (!modelChain %in% modelChain_vec ){print(paste0("No historical data for: ", modelChain));next}

		# read in data
		nc=read.csv(rcp26_file, sep=',', header=T, na.strings = "-999")
		
		# extract variable
		tas =nc[,modindex]
		if(all(is.na(tas)==TRUE)){print(paste0("no data found in variable, skipping this file", rcp26_file));next}

		greg_cal_rcp = strptime( nc[,1],format="%Y-%m-%d %H:%M:%S")

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
	save(df,file = paste0(outdir,"rcp26_",var))


	# Cordex rcp85 2006-2100=======================================================================
	rcp85_qmap_list <- list()
	rcp85_nqmap <- list()
	rcp85_qmap_season_list <- list()
	for (rcp85_file in rcp85_files){

		print(rcp85_file)

		modelNameBase = unlist(strsplit(rcp85_file,'/'))[length(unlist(strsplit(rcp85_file,'/')))]
		GCM = unlist(strsplit(modelNameBase,'_rcp85_'))[1]
		RCM = unlist(strsplit(modelNameBase,'_'))[4]
		modelChain=paste0(GCM,'_',RCM)

		# kip to next model if does not have hist data
		if (!modelChain %in% modelChain_vec ){print(paste0("No historical data for: ", modelChain));next}

		# define variable and gridbox
		nc=read.csv(rcp85_file, sep=',', header=T, na.strings = "-999")
		
		tas =nc[,modindex]
		if(all(is.na(tas)==TRUE)){print(paste0("no data found in variable, skipping this file", rcp85_file));next}

		greg_cal_rcp = strptime( nc[,1],format="%Y-%m-%d %H:%M:%S")

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
	df$Date<-greg_cal_rcp[si:ei]
	save(df,file = paste0(outdir,"rcp85_",var))



if(plot==TRUE){
#===============================================================================
# Plots
#===============================================================================
mycol = wes_palette("Zissou1", 3, type = "continuous")

#===============================================================================
# prepare data yearly
#===============================================================================

# extract data from list
hist_qmap_list_NOPARS = lapply(hist_qmap_season_list, `[[`, 1)
# prepare data annual means and stats
df=as.data.frame(hist_qmap_list_NOPARS)
years_hist=( substr(cordex_dates_cp, 1, 4)) 
hist_year = aggregate(df, list(years_hist), mean)
hist_year$SD=apply(hist_year[,2:dim(hist_year)[2]],1, sd, na.rm = TRUE)
hist_year$MEAN=apply(hist_year[,2:(dim(hist_year)[2]-1)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(hist_nqmap)
years_hist=( substr(cordex_dates_cp, 1, 4)) 
hist_year_nq = aggregate(df, list(years_hist), mean)
hist_year_nq$SD=apply(hist_year_nq[,2:dim(hist_year_nq)[2]],1, sd, na.rm = TRUE)
hist_year_nq$MEAN=apply(hist_year_nq[,2:(dim(hist_year_nq)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp26_qmap_season_list)
years_rcp26=( substr( greg_cal_rcp[si:ei], 1, 4)) 
rcp26_year = aggregate(df, list(years_rcp26), mean)
rcp26_year$SD=apply(rcp26_year[,2:dim(rcp26_year)[2]],1, sd, na.rm = TRUE)
rcp26_year$MEAN=apply(rcp26_year[,2:(dim(rcp26_year)[2]-1)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp26_nqmap)
years_rcp26=( substr( greg_cal_rcp[si:ei], 1, 4)) 
rcp26_year_nq = aggregate(df, list(years_rcp26), mean)
rcp26_year_nq$SD=apply(rcp26_year_nq[,2:dim(rcp26_year_nq)[2]],1, sd, na.rm = TRUE)
rcp26_year_nq$MEAN=apply(rcp26_year_nq[,2:(dim(rcp26_year_nq)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp85_qmap_season_list)
years_rcp85=( substr( greg_cal_rcp[si:ei], 1, 4)) 
rcp85_year = aggregate(df, list(years_rcp85), mean)
rcp85_year$SD=apply(rcp85_year[,2:dim(rcp85_year)[2]],1, sd, na.rm = TRUE)
rcp85_year$MEAN=apply(rcp85_year[,2:(dim(rcp85_year)[2]-1)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp85_nqmap)
years_rcp85=( substr( greg_cal_rcp[si:ei], 1, 4)) 
rcp85_year_nq = aggregate(df, list(years_rcp85), mean)
rcp85_year_nq$SD=apply(rcp85_year_nq[,2:dim(rcp85_year_nq)[2]],1, sd, na.rm = TRUE)
rcp85_year_nq$MEAN=apply(rcp85_year_nq[,2:(dim(rcp85_year_nq)[2]-1)],1, mean, na.rm = TRUE)

#obs
years_obs=( substr( obs_cp_dates, 1, 4)) 
obs_year = aggregate(obs_cp_period, list(years_obs), mean)

#===envelope plots========================================================================
if (var == 'TA'){

	pdf(paste(outdir,var,"TS2.pdf"))
	# caption: Mean annual near surface air temperature at the Weissfluhjoch (2540m asl)  showing corrected historical, RCP2.6 and RCP8.5 timseries. Observations and uncorrected historical data are also shown for comparison. The coloured envelops indicate +/- 1 SD of the model spread and multi-modal mean is given by the bold line.
	lwd=3
	plot(rcp26_year$Group.1, rcp26_year$MEAN, xlim=c(1979,2100),ylim=c(270,280), type='l', col=mycol[1], lwd=lwd, ylab="Air temperacture (K)", xlab=" ")


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

	lines(hist_year_nq, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd,lty=2)
	lines(obs_year,lwd=lwd)
	legend("bottomright", c("Hist", "RCP2.6", "RCP8.5", "Hist_uncorrected", "OBS"), col=c(rgb(0, 1, 0),mycol[1],mycol[3],mycol[2], "black"),lty=c(1,1,1,2,1),lwd=3)
	dev.off()

}

if (var == 'TA'){

	pdf(paste(outdir,var,"TS2_NOQMAP.pdf"))
	# caption: Mean annual near surface air temperature at the Weissfluhjoch (2540m asl)  showing corrected historical, RCP2.6 and RCP8.5 timseries. Observations and uncorrected historical data are also shown for comparison. The coloured envelops indicate +/- 1 SD of the model spread and multi-modal mean is given by the bold line.
	lwd=3
	plot(rcp26_year$Group.1, rcp26_year$MEAN, xlim=c(1979,2100),ylim=c(270,280), type='l', col=mycol[1], lwd=lwd, ylab="Air temperacture (K)", xlab=" ")
	lines(rcp26_year_nq$Group.1, rcp26_year_nq$MEAN, ylim=c(270,280), type='l', col=mycol[1],lwd=lwd,lty=2)

	lines(rcp85_year$Group.1, rcp85_year$MEAN, ylim=c(270,280), type='l', col=mycol[3],lwd=lwd)
	lines(rcp85_year_nq$Group.1, rcp85_year_nq$MEAN, ylim=c(270,280), type='l', col=mycol[3],lwd=lwd,lty=2)


	lines(hist_year$Group.1, hist_year$MEAN, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd)
	lines(hist_year_nq, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd,lty=2)
	lines(obs_year,lwd=lwd)
	legend("bottomright", c("Hist", "RCP2.6", "RCP8.5", "Hist_uncorrected", "OBS"), col=c(rgb(0, 1, 0),mycol[1],mycol[3],mycol[2], "black"),lty=c(1,1,1,2,1),lwd=3)
	dev.off()

}


#===============================================================================
# prepare data monthly
#===============================================================================

# extract data from list
hist_qmap_season_list_NOPARS = lapply(hist_qmap_season_list, `[[`, 1)
# prepare data annual means and stats
df=as.data.frame(hist_qmap_season_list_NOPARS)
months_hist=( substr(cordex_dates_cp, 1, 7)) 
hist_month = aggregate(df, list(months_hist), mean)
hist_month$SD=apply(hist_month[,2:dim(hist_month)[2]],1, sd, na.rm = TRUE)
hist_month$MEAN=apply(hist_month[,2:(dim(hist_month)[2]-1)],1, mean, na.rm = TRUE)


# extract data from list
hist_qmap_list_NOPARS = lapply(hist_qmap_list, `[[`, 1)
# no season qmap
df=as.data.frame(hist_qmap_list_NOPARS)
months_hist=( substr(cordex_dates_cp, 1, 7)) 
hist_month_noseas = aggregate(df, list(months_hist), mean)
hist_month_noseas $SD=apply(hist_month_noseas [,2:dim(hist_month_noseas )[2]],1, sd, na.rm = TRUE)
hist_month_noseas $MEAN=apply(hist_month_noseas [,2:(dim(hist_month_noseas )[2]-1)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(hist_nqmap)
months_hist=( substr(cordex_dates_cp, 1, 7)) 
hist_month_nq = aggregate(df, list(months_hist), mean)
hist_month_nq$SD=apply(hist_month_nq[,2:dim(hist_month_nq)[2]],1, sd, na.rm = TRUE)
hist_month_nq$MEAN=apply(hist_month_nq[,2:(dim(hist_month_nq)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp26_qmap_season_list)
months_rcp26=( substr( greg_cal_rcp[si:ei], 1, 7)) 
rcp26_month = aggregate(df, list(months_rcp26), mean)
rcp26_month$SD=apply(rcp26_month[,2:dim(rcp26_month)[2]],1, sd, na.rm = TRUE)
rcp26_month$MEAN=apply(rcp26_month[,2:(dim(rcp26_month)[2]-1)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp26_nqmap)
months_rcp26=( substr( greg_cal_rcp[si:ei], 1, 7)) 
rcp26_month_nq = aggregate(df, list(months_rcp26), mean)
rcp26_month_nq$SD=apply(rcp26_month_nq[,2:dim(rcp26_month_nq)[2]],1, sd, na.rm = TRUE)
rcp26_month_nq$MEAN=apply(rcp26_month_nq[,2:(dim(rcp26_month_nq)[2]-1)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp85_qmap_season_list)
months_rcp85=( substr( greg_cal_rcp[si:ei], 1, 7)) 
rcp85_month = aggregate(df, list(months_rcp85), mean)
rcp85_month$SD=apply(rcp85_month[,2:dim(rcp85_month)[2]],1, sd, na.rm = TRUE)
rcp85_month$MEAN=apply(rcp85_month[,2:(dim(rcp85_month)[2]-1)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp85_nqmap)
months_rcp85=( substr( greg_cal_rcp[si:ei], 1, 7)) 
rcp85_month_nq = aggregate(df, list(months_rcp85), mean)
rcp85_month_nq$SD=apply(rcp85_month_nq[,2:dim(rcp85_month_nq)[2]],1, sd, na.rm = TRUE)
rcp85_month_nq$MEAN=apply(rcp85_month_nq[,2:(dim(rcp85_month_nq)[2]-1)],1, mean, na.rm = TRUE)

#obs
months_obs=( substr( obs_cp_dates, 1, 7)) 
obs_month = aggregate(obs_cp_period, list(months_obs), mean)



#===============================================================================
# summary stats
#===============================================================================
# cor(obs_month, hist_month$MEAN)
# cor(obs_month, hist_qmap_month$MEAN)

# rmse(obs_month, hist_month$MEAN)
# rmse(obs_month, hist_qmap_month$MEAN)

# mean(obs_month)
# mean(hist_month$MEAN)
# mean(hist_qmap_month$MEAN)

#===============================================================================
# ecdf
#===============================================================================
lwd=3
plot(ecdf(obs_month$x), ylab="cdf",  col=mycol[1], xlab="TA (K)", main=var, lwd=lwd,lty=1)
lines(ecdf(hist_month_nq[,2]),  col=mycol[2], lwd=lwd)
lines(ecdf(hist_month_noseas[,2]),  col=mycol[3], lwd=1)
lines(ecdf(hist_month[,2]),  col='green', lwd=1)
legend("topright", c("OBS", "CLIM", "CLIM_QM", "CLIM_QMSEASON"), col=c(mycol[1],mycol[2], mycol[3], 'green'),lty=c(1,1,1,1) ,lwd=lwd)


#===============================================================================
# doy
#===============================================================================

# make doy vector (slightly different periods (+/- 1 year))
doy_hist=substr((cordex_dates_cp), 6,10)
doy_obs=substr((obs_cp_dates), 6,10)

# DOY obs
obs_doy = aggregate(obs_cp_period, list(doy_obs), mean)

# DOY SEASON QMAP
hist_qmap_season_list_NOPARS = lapply(hist_qmap_season_list, `[[`, 1)
df=as.data.frame(hist_qmap_season_list_NOPARS)
hist_doy = aggregate(df, list(doy_hist), mean)
hist_doy$MEAN=apply(hist_doy[,2:(dim(hist_doy)[2]-1)],1, mean, na.rm = TRUE)

# DOY NO SEASON QMAP
hist_qmap_list_NOPARS = lapply(hist_qmap_list, `[[`, 1)
df=as.data.frame(hist_qmap_list_NOPARS)
hist_doy_noseas = aggregate(df, list(doy_hist), mean)
hist_doy_noseas $MEAN=apply(hist_doy_noseas [,2:(dim(hist_doy_noseas )[2]-1)],1, mean, na.rm = TRUE)

# NO QMAP
df=as.data.frame(hist_nqmap)
hist_doy_nq = aggregate(df, list(doy_hist), mean)
hist_doy_nq$MEAN=apply(hist_doy_nq[,2:(dim(hist_doy_nq)[2]-1)],1, mean, na.rm = TRUE)


#
#pdf(paste(outdir,var,"DOY2.pdf"))
plot(obs_doy$x, type='l',col=mycol[1],lwd=lwd,lty=2, xlab='DOY', ylab='Air temperature [k]', main=var)
lines(hist_doy_nq$MEAN, type='l', col=mycol[2], lwd=lwd)
lines(hist_doy_noseas$MEAN ,type='l', col=mycol[3],lwd=lwd)
lines(hist_doy$MEAN, type='l', col='green',lwd=lwd)
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



}

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

outdir=paste0(indir,"/fsm/")
results = list.files(path = "/home/joel/sim/qmap/topoclim_test/aqmap_results/",full.names=T )
hist = results[grep("hist",results)]
rcp26 = results[grep("rcp26",results)]
rcp85 = results[grep("rcp85",results)]

#DW

typ='HIST'
load(hist[1] )
DW<-df
load(hist[2] )
ILWR<-df
load(hist[3] )
ISWR<-df
load(hist[4] )
P<-df
load(hist[5] )
PINT<-df
load(hist[6] )
RH<-df
load(hist[7] )
TA<-df
load(hist[8] )
VW<-df

models = names(df)[1:length(names(df))-1]
dates = df[length(df)]
for (model in models){
	print(model)
	modellist=list()
	modellist['datetime'] <- dates
	modellist['DW']<- DW[model]
	modellist['ILWR']<- ILWR[model]
	modellist['ISWR']<- ISWR[model] # qmap can make non-zero night time rad which is unphysical
	modellist['P']<- P[model]
	modellist['PINT']<- PINT[model]
	modellist['RH']<- RH[model]
	modellist['TA']<- TA[model]
	modellist['VW']<- VW[model]
	df = as.data.frame(modellist)
	df[c(2:5, 7:9)] =round(df[c(2:5, 7:9)],1)
	df[6 ] = round(df[c(6)],4)
	write.csv(df, paste0(outdir, model, "_",typ,"_QMAP.txt"), row.names=FALSE)
	
}


typ='RCP26'
load(rcp26[1] )
DW<-df
load(rcp26[2] )
ILWR<-df
load(rcp26[3] )
ISWR<-df
load(rcp26[4] )
P<-df
load(rcp26[5] )
PINT<-df
load(rcp26[6] )
RH<-df
load(rcp26[7] )
TA<-df
load(rcp26[8] )
VW<-df

models = names(df)[1:length(names(df))-1]
dates = df[length(df)]
for (model in models){
	print(model)
	modellist=list()
	modellist['datetime'] <- dates
	modellist['DW']<- DW[model]
	modellist['ILWR']<- ILWR[model]
	modellist['ISWR']<- ISWR[model]
	modellist['P']<- P[model]
	modellist['PINT']<- PINT[model]
	modellist['RH']<- RH[model]
	modellist['TA']<- TA[model]
	modellist['VW']<- VW[model]
	df = as.data.frame(modellist)
	df[c(2:5, 7:9)] =round(df[c(2:5, 7:9)],1)
	df[6 ] = round(df[c(6)],4)
	write.csv(df, paste0(outdir, model, "_",typ,"_QMAP.txt"), row.names=FALSE)
	
}


typ='RCP85'
load(rcp85[1] )
DW<-df
load(rcp85[2] )
ILWR<-df
load(rcp85[3] )
ISWR<-df
load(rcp85[4] )
P<-df
load(rcp85[5] )
PINT<-df
load(rcp85[6] )
RH<-df
load(rcp85[7] )
TA<-df
load(rcp85[8] )
VW<-df



models = names(df)[1:length(names(df))-1]
dates = df[length(df)]
for (model in models){
	print(model)
	modellist=list()
	modellist['datetime'] <- dates
	modellist['DW']<- DW[model]
	modellist['ILWR']<- ILWR[model]
	modellist['ISWR']<- ISWR[model] 
	modellist['P']<- P[model]
	modellist['PINT']<- PINT[model]
	modellist['RH']<- RH[model]
	modellist['TA']<- TA[model]
	modellist['VW']<- VW[model]
	df = as.data.frame(modellist)
	df[c(2:5, 7:9)] =round(df[c(2:5, 7:9)],1)
	df[6 ] = round(df[c(6)],4)
	write.csv(df, paste0(outdir, model, "_",typ,"_QMAP.txt"), row.names=FALSE)
	
}

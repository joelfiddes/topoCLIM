args = commandArgs(trailingOnly=TRUE)
wd = args[1] # wd = '/home/joel/sim/qmap/tclim_points_paper/'
sample = args[2]#189
daily_obs = args[3]


# meteo
monthly_met = "../examples/eval/wfj_station_monthly.csv"
doy_met = "../examples/eval/wfj_station_doy.csv"
monthly_met = read.csv(monthly_met)
doy_met = read.csv(doy_met)


#indir = "/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc1_1D/"
#daily_obs  = '/home/joel/sim/qmap/GR_data/sim/g4//forcing/meteoc1_1D.csv' 

indir=paste0(wd, '/s',sample,'/')

pdf(paste0(indir,"evalplot_singleModel_season_REV.pdf"),height=20, width=10)
par(mfrow=c(5,2))
par(mar=c(5,5,4,2))
require(wesanderson)




# the period to do quantile mapping over (obs and sim must overlap!)
startDateQmap = "1990-01-01"
endDateQmap = "1999-12-31"

# common historical period of models
startDateHist = "1980-01-01"
endDateHist = "2005-12-31"


startDateQmap = "1990-01-01"
endDateQmap = "1999-12-31"
# common historical period of models
startDateHist = "1980-01-01"
endDateHist = "2005-12-31"

# obs

obs=read.csv(daily_obs)
	# aggregate obs to daily resolution

qmap_files = list.files(path = paste0(indir,"/aqmap_results/"),full.names=TRUE )
results_qmap = qmap_files[ grep('_qmap_',qmap_files)]
results_nmap = qmap_files[ grep('_nmap_',qmap_files)]
results_season = qmap_files[ grep('_season_',qmap_files)]

# stats vectors
rMat=c()
rmseMat=c()
biasMat=c()

rMat1=c()
rmseMat1=c()
biasMat1=c()

# llop thru list of vars
myvars=c('tas' ,'pr', 'rsds','rlds' , 'hurs')#, 'ps','uas')

for (var in myvars){
print(var)


	# if(var=="tas"){obsindex <-11; convFact <-1} # K
	# if(var=="tasmin"){obsindex <-14; convFact <-1} # K
	# if(var=="tasmax"){obsindex <-13; convFact <-1} # K
	# 	if(var=="pr"){obsindex <-6; convFact <-(1/3600)} # kgm-2s-1 obs are in mean mm/hr for that day convert to mm/s (cordex)
	# 		if(var=="ps"){obsindex <-5; convFact <-1}# Pa
	# 			if(var=="hurs"){obsindex <-8; convFact <-100} # % 0-100
	# 				if(var=="rsds"){obsindex <-4; convFact <-1}	# Wm-2
	# 					if(var=="rlds"){obsindex <-3; convFact <-1}# Wm-2
	# 						if(var=="uas"){obsindex <-2; convFact <-1}# ms-1
	# 							if(var=="vas"){obsindex <-2; convFact <-1}# ms-1

	if(var=="tas"){obsindex <-6; convFact <-1; xlab <- "Mean daily air temperature [K]"} # K
	if(var=="tasmin"){obsindex <-11; convFact <-1;  xlab <- "Min daily air temperature [K]"} # K
	if(var=="tasmax"){obsindex <-10; convFact <-1 ;  xlab <- "Max daily air temperature [K]"} # K
		if(var=="pr"){obsindex <-12; convFact <-(1);  xlab <- "Mean daily precipitation rate [Kg m^-2 s^-1]"} # kgm-2s-1 obs are in mean mm/hr for that day convert to mm/s (cordex)
			if(var=="ps"){obsindex <-9; convFact <-1 ;  xlab <- "Mean daily air pressure [Pa]" }# Pa
				if(var=="hurs"){obsindex <-7; convFact <-1;  xlab <- "Mean daily relative humidity [%]" } # % 0-100
					if(var=="rsds"){obsindex <-2; convFact <-1 ;  xlab <- "Mean daily incoming shortwave radiation [W m-2]"}	# Wm-2
						if(var=="rlds"){obsindex <-3; convFact <-1 ;  xlab <- "Mean daily incoming longwave radiation [W m-2]"}# Wm-2
							if(var=="uas"){obsindex <-8; convFact <-1}# ms-1
								if(var=="vas"){obsindex <-8; convFact <-1}# ms-1
								# where is DW?

 	if(var=="tas"){obsindex <-6;metindex <-2; convFact <-1;convFactMet<- 1; xlab <- "Mean daily air temperature [K]"} # K
	if(var=="tasmin"){obsindex <-11; convFact <-1;convFactMet<- 1;  xlab <- "Min daily air temperature [K]"} # K
	if(var=="tasmax"){obsindex <-10; convFact <-1;convFactMet<- 1 ;  xlab <- "Max daily air temperature [K]"} # K
		if(var=="pr"){obsindex <-12;metindex <-10; convFact <-(1);convFactMet<- 0.000555; xlab <-  ~paste("Mean daily precipitation rate [Kg m"^-2 ,"s"^-1,"]")} # expression(Speed ~ ms^-1 ~ by ~ impeller)    kgm-2s-1 obs are in mean mm/hr for that day convert to mm/s (cordex)
			
			if(var=="ps"){obsindex <-9; metindex <-3;convFact <-1;convFactMet<- 1 ;  xlab <- "Mean daily air pressure [Pa]" }# Pa
				if(var=="hurs"){obsindex <-7;metindex <-3 ;convFact <-1;convFactMet<- 100;  xlab <- "Mean daily relative humidity [%]" } # % 0-100
					if(var=="rsds"){obsindex <-2; metindex <-6;convFact <-1;convFactMet<- 1 ;  xlab <- ~paste("Mean daily incoming shortwave radiation [W m"^-2,"]")}	# Wm-2
						if(var=="rlds"){obsindex <-3; metindex <-8;convFact <-1;convFactMet<- 1 ;  xlab <- ~paste("Mean daily incoming longwave radiation [W m"^-2,"]")}# Wm-2
							if(var=="uas"){obsindex <-8; convFact <-1;convFactMet<- 1}# ms-1
								if(var=="vas"){obsindex <-8; convFact <-1;convFactMet<- 1}# ms-1
								# where is DW?

   #read met
  doy_metPar= doy_met[,metindex]*convFactMet
  month_metPar= monthly_met[,metindex]*convFactMet

require(caTools)
  if (var =="pr"){

  	doy_metPar = runmean(doy_metPar,10)
  }

	# read and convert obs unit
	obs_var=obs[,obsindex]*convFact	

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

# filenames
load(results_season [grep(paste0('hist.*',var,'$'),results_season) ] )
hist_qmap_season_list = df
load(results_nmap [grep(paste0('hist.*',var,'$'),results_season) ] )
hist_nqmap= df
load(results_qmap [grep(paste0('hist.*',var,'$'),results_season) ] )
hist_qmap_list= df

load(results_season [grep(paste0('rcp26.*',var,'$'),results_season) ] )
rcp26_qmap_season_list= df
load(results_nmap [grep(paste0('rcp26.*',var,'$'),results_season) ] )
rcp26_nqmap=df
load(results_season [grep(paste0('rcp85.*',var,'$'),results_season) ] )
rcp85_qmap_season_list= df
load(results_nmap [grep(paste0('rcp85.*',var,'$'),results_season) ] )
rcp85_nqmap=df

# do stats here already


cordex_dates_cp=hist_qmap_season_list$Date
greg_cal_rcp=rcp26_qmap_season_list$Date
#===============================================================================
# Plots
#===============================================================================
mycol = wes_palette("Zissou1", 3, type = "continuous")

#===============================================================================
# prepare data yearly
#===============================================================================

# extract data from list

# prepare data annual means and stats
df=as.data.frame(hist_qmap_season_list)
years_hist=( substr(cordex_dates_cp, 1, 4)) 
hist_year = aggregate(df, list(years_hist), mean)
hist_year$SD=apply(hist_year[,2:(dim(hist_year)[2]-1)],1, sd, na.rm = TRUE)
hist_year$MEAN=apply(hist_year[,2:(dim(hist_year)[2]-2)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(hist_nqmap)
years_hist=( substr(cordex_dates_cp, 1, 4)) 
hist_year_nq = aggregate(df, list(years_hist), mean)
hist_year_nq$SD=apply(hist_year_nq[,2:(dim(hist_year_nq)[2]-1)],1, sd, na.rm = TRUE)
hist_year_nq$MEAN=apply(hist_year_nq[,2:(dim(hist_year_nq)[2]-2)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp26_qmap_season_list)
years_rcp26=( substr( greg_cal_rcp, 1, 4)) 
rcp26_year = aggregate(df, list(years_rcp26), mean)
rcp26_year$SD=apply(rcp26_year[,2:(dim(rcp26_year)[2]-1)],1, sd, na.rm = TRUE)
rcp26_year$MEAN=apply(rcp26_year[,2:(dim(rcp26_year)[2]-2)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp26_nqmap)
years_rcp26=( substr( greg_cal_rcp, 1, 4)) 
rcp26_year_nq = aggregate(df, list(years_rcp26), mean)
rcp26_year_nq$SD=apply(rcp26_year_nq[,2:(dim(rcp26_year_nq)[2]-1)],1, sd, na.rm = TRUE)
rcp26_year_nq$MEAN=apply(rcp26_year_nq[,2:(dim(rcp26_year_nq)[2]-2)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp85_qmap_season_list)
years_rcp85=( substr( greg_cal_rcp, 1, 4)) 
rcp85_year = aggregate(df, list(years_rcp85), mean)
rcp85_year$SD=apply(rcp85_year[,2:(dim(rcp85_year)[2]-1)],1, sd, na.rm = TRUE)
rcp85_year$MEAN=apply(rcp85_year[,2:(dim(rcp85_year)[2]-2)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp85_nqmap)
years_rcp85=( substr( greg_cal_rcp, 1, 4)) 
rcp85_year_nq = aggregate(df, list(years_rcp85), mean)
rcp85_year_nq$SD=apply(rcp85_year_nq[,2:(dim(rcp85_year_nq)[2]-1)],1, sd, na.rm = TRUE)
rcp85_year_nq$MEAN=apply(rcp85_year_nq[,2:(dim(rcp85_year_nq)[2]-2)],1, mean, na.rm = TRUE)

#obs
years_obs=( substr( obs_cp_dates, 1, 4)) 
obs_year = aggregate(obs_cp_period, list(years_obs), mean)

#===envelope plots========================================================================
if (var == 'tas'){

	pdf(paste0(indir,var,"TS2_REV.pdf"))
	# caption: Mean annual near surface air temperature at the Weissfluhjoch (2540m asl)  showing corrected historical, RCP2.6 and RCP8.5 timseries. Observations and uncorrected historical data are also shown for comparison. The coloured envelops indicate +/- 1 SD of the model spread and multi-modal mean is given by the bold line.
	lwd=3
	plot(rcp26_year$Group.1, rcp26_year$MEAN-273.15, xlim=c(1979,2100),ylim=c(270-273.15,280-273.15), type='l', col="white", lwd=lwd, ylab="Air temperature (°C)", xlab=" ")


 	y = c(rcp26_year$MEAN-273.15- rcp26_year$SD,rev(rcp26_year$MEAN-273.15+ rcp26_year$SD))
 	x=c((rcp26_year$Group.1),rev(rcp26_year$Group.1))
 	polygon (x,y, col=rgb(0, 0, 1,0.3),border='NA')
 lines(rcp26_year$Group.1, rcp26_year$MEAN-273.15, ylim=c(270-273.15,280-273.15), type='l', col=mycol[1],lwd=lwd,lty=1)
 lines(rcp26_year_nq$Group.1, rcp26_year_nq$MEAN-273.15, ylim=c(270-273.15,280-273.15), type='l', col=mycol[1],lwd=lwd,lty=2)

 	lines(rcp85_year$Group.1, rcp85_year$MEAN-273.15, ylim=c(270-273.15,280-273.15), type='l', col=mycol[3],lwd=lwd)
 	y = c(rcp85_year$MEAN-273.15- rcp85_year$SD,rev(rcp85_year$MEAN-273.15+ rcp85_year$SD))
 	x=c((rcp85_year$Group.1),rev(rcp85_year$Group.1))
 	polygon (x,y, col=rgb(1, 0, 0,0.3),border='NA')
 lines(rcp85_year_nq$Group.1, rcp85_year_nq$MEAN-273.15, ylim=c(270-273.15,280-273.15), type='l', col=mycol[3],lwd=lwd,lty=2)

	
	y = c(hist_year$MEAN-273.15- hist_year$SD,rev(hist_year$MEAN-273.15+ hist_year$SD))
 	x=c((hist_year$Group.1),rev(hist_year$Group.1))
 	polygon (x,y, col=rgb(0, 1, 0,0.3),border='NA')
 	lines(hist_year$Group.1, hist_year$MEAN-273.15, ylim=c(270-273.15,280-273.15), type='l', col=mycol[2],lwd=lwd)


	lines(hist_year$Group.1, hist_year_nq$MEAN-273.15, type='l', col=mycol[2],lwd=lwd,lty=2)
 	lines(obs_year$Group.1 , obs_year$x-273.15,lwd=lwd, lty=1, col="grey31")
 	abline(h=0, col="grey31")
 	legend("topleft", c("Hist","Hist_raw", "RCP2.6", "RCP2.6_raw", "RCP8.5" ,"RCP8.5_raw",  "T-MET"), col=c(mycol[2],mycol[2], mycol[1],mycol[1],mycol[3],mycol[3], "grey31"),lty=c(1,2,1,2,1,2,1),lwd=3)
	dev.off()

}

if (var == 'tas'){

	pdf(paste0(indir,var,"TS2_NOQMAP_REV.pdf"))
	# caption: Mean annual near surface air temperature at the Weissfluhjoch (2540m asl)  showing corrected historical, RCP2.6 and RCP8.5 timseries. Observations and uncorrected historical data are also shown for comparison. The coloured envelops indicate +/- 1 SD of the model spread and multi-modal mean is given by the bold line.
	lwd=3
	plot(rcp26_year$Group.1, rcp26_year$MEAN, xlim=c(1979,2100),ylim=c(270,280), type='l', col=mycol[1], lwd=lwd, ylab="Air temperature (K)", xlab=" ")
	lines(rcp26_year_nq$Group.1, rcp26_year_nq$MEAN, ylim=c(270,280), type='l', col=mycol[1],lwd=lwd,lty=2)

	lines(rcp85_year$Group.1, rcp85_year$MEAN, ylim=c(270,280), type='l', col=mycol[3],lwd=lwd)
	lines(rcp85_year_nq$Group.1, rcp85_year_nq$MEAN, ylim=c(270,280), type='l', col=mycol[3],lwd=lwd,lty=2)


	lines(hist_year$Group.1, hist_year$MEAN, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd)
	lines(hist_year_nq, ylim=c(270,280), type='l', col=mycol[2],lwd=lwd,lty=2)
	lines(obs_year,lwd=lwd)
	legend("bottomright", c("Hist", "RCP2.6", "RCP8.5", "Hist_uncorrected", "T-MET"), col=c(rgb(0, 1, 0),mycol[1],mycol[3],mycol[2], "black"),lty=c(1,1,1,2,1),lwd=3)
	dev.off()

}


#===============================================================================
# prepare data monthly
#===============================================================================

# extract data from list

# prepare data annual means and stats
df=as.data.frame(hist_qmap_season_list)
months_hist=( substr(cordex_dates_cp, 1, 7)) 
hist_month = aggregate(df, list(months_hist), mean)
hist_month$SD=apply(hist_month[,2:(dim(hist_month)[2]-1)],1, sd, na.rm = TRUE)
hist_month$MEAN=apply(hist_month[,2:(dim(hist_month)[2]-2)],1, mean, na.rm = TRUE)


# extract data from list

# no season qmap
df=as.data.frame(hist_qmap_list)
months_hist=( substr(cordex_dates_cp, 1, 7)) 
hist_month_noseas = aggregate(df, list(months_hist), mean)
hist_month_noseas $SD=apply(hist_month_noseas [,2:(dim(hist_month_noseas )[2]-1)],1, sd, na.rm = TRUE)
hist_month_noseas $MEAN=apply(hist_month_noseas [,2:(dim(hist_month_noseas )[2]-2)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(hist_nqmap)
months_hist=( substr(cordex_dates_cp, 1, 7)) 
hist_month_nq = aggregate(df, list(months_hist), mean)
hist_month_nq$SD=apply(hist_month_nq[,2:(dim(hist_month_nq)[2]-1)],1, sd, na.rm = TRUE)
hist_month_nq$MEAN=apply(hist_month_nq[,2:(dim(hist_month_nq)[2]-2)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp26_qmap_season_list)
months_rcp26=( substr( greg_cal_rcp, 1, 7)) 
rcp26_month = aggregate(df, list(months_rcp26), mean)
rcp26_month$SD=apply(rcp26_month[,2:(dim(rcp26_month)[2]-1)],1, sd, na.rm = TRUE)
rcp26_month$MEAN=apply(rcp26_month[,2:(dim(rcp26_month)[2]-2)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp26_nqmap)
months_rcp26=( substr( greg_cal_rcp, 1, 7)) 
rcp26_month_nq = aggregate(df, list(months_rcp26), mean)
rcp26_month_nq$SD=apply(rcp26_month_nq[,2:(dim(rcp26_month_nq)[2]-1)],1, sd, na.rm = TRUE)
rcp26_month_nq$MEAN=apply(rcp26_month_nq[,2:(dim(rcp26_month_nq)[2]-2)],1, mean, na.rm = TRUE)

df=as.data.frame(rcp85_qmap_season_list)
months_rcp85=( substr( greg_cal_rcp, 1, 7)) 
rcp85_month = aggregate(df, list(months_rcp85), mean)
rcp85_month$SD=apply(rcp85_month[,2:(dim(rcp85_month)[2]-1)],1, sd, na.rm = TRUE)
rcp85_month$MEAN=apply(rcp85_month[,2:(dim(rcp85_month)[2]-2)],1, mean, na.rm = TRUE)

# no qmap
df=as.data.frame(rcp85_nqmap)
months_rcp85=( substr( greg_cal_rcp, 1, 7)) 
rcp85_month_nq = aggregate(df, list(months_rcp85), mean)
rcp85_month_nq$SD=apply(rcp85_month_nq[,2:(dim(rcp85_month_nq)[2]-1)],1, sd, na.rm = TRUE)
rcp85_month_nq$MEAN=apply(rcp85_month_nq[,2:(dim(rcp85_month_nq)[2]-2)],1, mean, na.rm = TRUE)

#obs
months_obs=( substr( obs_cp_dates, 1, 7)) 
obs_month = aggregate(obs_cp_period, list(months_obs), mean)






#===============================================================================
# ecdf
#===============================================================================
lwd=2
plot(ecdf(obs_month$x), ylab="cdf",  col="black", xlab=xlab, main=var, lwd=lwd,lty=1,  cex.main=2, cex.lab=1.2)
lines(ecdf(hist_month_nq[,2]),  col=mycol[2], lwd=lwd)
lines(ecdf(hist_month_noseas[,2]),  col=mycol[1], lwd=lwd)
lines(ecdf(hist_month[,2]), col="orange", lwd=lwd)
lines(ecdf(month_metPar), col="darkgrey")
legend("topleft", c("STATION","T-MET", "CLIM", "QM", "QM_MONTH"), col=c("darkgrey","black", mycol[2], mycol[1], "orange"),lty=c(1,1,1,1,1) ,lwd=lwd,bg = "white")


#===============================================================================
# doy
#===============================================================================

# make doy vector (slightly different periods (+/- 1 year))
doy_hist=substr((cordex_dates_cp), 6,10)
doy_obs=substr((obs_cp_dates), 6,10)

# DOY obs
obs_doy = aggregate(obs_cp_period, list(doy_obs), mean)

# DOY SEASON QMAP

df=as.data.frame(hist_qmap_season_list)
hist_doy = aggregate(df, list(doy_hist), mean)
hist_doy$MEAN=apply(hist_doy[,2:(dim(hist_doy)[2]-1)],1, mean, na.rm = TRUE)

# DOY NO SEASON QMAP

df=as.data.frame(hist_qmap_list)
hist_doy_noseas = aggregate(df, list(doy_hist), mean)
hist_doy_noseas $MEAN=apply(hist_doy_noseas [,2:(dim(hist_doy_noseas )[2]-1)],1, mean, na.rm = TRUE)

# NO QMAP
df=as.data.frame(hist_nqmap)
hist_doy_nq = aggregate(df, list(doy_hist), mean)
hist_doy_nq$MEAN=apply(hist_doy_nq[,2:(dim(hist_doy_nq)[2]-1)],1, mean, na.rm = TRUE)


#
#pdf(paste(outdir,var,"DOY2.pdf"))
ylim= c(min(c(obs_doy$x,hist_doy_nq$MEAN,hist_doy_noseas$MEAN,hist_doy$MEAN)) , max(c(obs_doy$x,hist_doy_nq$MEAN,hist_doy_noseas$MEAN,hist_doy$MEAN)))
plot(obs_doy$x, type='l',col="black",lwd=lwd,lty=2, xlab='DOY', ylab=xlab, main=var, ylim=ylim, cex.main=2, cex.lab=1.2)
lines(hist_doy_nq$MEAN, type='l', col=mycol[2], lwd=lwd)
lines(hist_doy_noseas$MEAN ,type='l', col=mycol[1],lwd=lwd)
lines(hist_doy$MEAN, type='l', col="orange",lwd=lwd)
lines(doy_metPar, col="darkgrey", lwd=2)
#legend("topright", c("OBS", "CLIM", "CLIM_QM"), col=c(mycol[1],mycol[2], mycol[3]),lty=c(2,1,1) ,lwd=lwd)
legend("topright", c("STATION","T-MET", "CLIM", "QM", "QM_MONTH"), col=c("darkgrey","black", mycol[2], mycol[1],"orange"),lty=c(1,2,1,1,1) ,lwd=lwd, bg = "white")



#===============================================================================
# summary stats
#===============================================================================
require(hydroGOF)
clim_r = cor(obs_doy$x, hist_doy_nq$MEAN)
qm_r = cor(obs_doy$x, hist_doy_noseas$MEAN)
qms_r = cor(obs_doy$x, hist_doy$MEAN)


clim_rm = rmse(hist_doy_nq$MEAN, obs_doy$x)
qm_rm = rmse(hist_doy_noseas$MEAN, obs_doy$x)
qms_rm = rmse(hist_doy$MEAN,  obs_doy$x)

clim_pb = pbias(hist_doy_nq$MEAN, obs_doy$x)
qm_pb = pbias(hist_doy_noseas$MEAN, obs_doy$x)
qms_pb = pbias(hist_doy$MEAN,  obs_doy$x)



rvec=c(clim_r, qm_r, qms_r)
rmsevec=c(clim_rm, qm_rm, qms_rm)
biasvec=c(clim_pb, qm_pb, qms_pb)

rMat=cbind(rMat,rvec)
rmseMat=cbind(rmseMat,rmsevec)
biasMat=cbind(biasMat, biasvec)

# against meteo
clim_r = cor(doy_metPar, hist_doy_nq$MEAN)
qm_r = cor(doy_metPar, hist_doy_noseas$MEAN)
qms_r = cor(doy_metPar, hist_doy$MEAN)
tmet_r = cor(doy_metPar, obs_doy$x)

clim_rm = rmse(hist_doy_nq$MEAN, doy_metPar)
qm_rm = rmse(hist_doy_noseas$MEAN, doy_metPar)
qms_rm = rmse(hist_doy$MEAN,  doy_metPar)
tmet_rm = rmse(obs_doy$x,doy_metPar)

clim_pb = pbias(hist_doy_nq$MEAN, doy_metPar)
qm_pb = pbias(hist_doy_noseas$MEAN, doy_metPar)
qms_pb = pbias(hist_doy$MEAN,  doy_metPar)
tmet_pb = pbias(obs_doy$x,doy_metPar)


rvec=c(clim_r, qm_r, qms_r, tmet_r)
rmsevec=c(clim_rm, qm_rm, qms_rm, tmet_rm)
biasvec=c(clim_pb, qm_pb, qms_pb, tmet_pb)

rMat1=cbind(rMat1,rvec)
rmseMat1=cbind(rmseMat1,rmsevec)
biasMat1=cbind(biasMat1, biasvec)

}
dev.off()

df = data.frame(rMat)
colnames(df) <- c('tas' ,'pr', 'rsds','rlds' , 'hurs')
write.table(df, paste0(wd, "/r_stats.csv"))

df = data.frame(rmseMat)
colnames(df) <- c('tas' ,'pr', 'rsds','rlds' , 'hurs')
write.table(df, paste0(wd, "/rmse_stats.csv"))

df = data.frame(biasMat)
colnames(df) <- c('tas' ,'pr', 'rsds','rlds' , 'hurs')
write.table(df, paste0(wd, "/bias_stats.csv"))

df = data.frame(rMat1)
colnames(df) <- c('tas' ,'pr', 'rsds','rlds' , 'hurs')
write.table(df, paste0(wd, "/r_stats1.csv"))

df = data.frame(rmseMat1)
colnames(df) <- c('tas' ,'pr', 'rsds','rlds' , 'hurs')
write.table(df, paste0(wd, "/rmse_stats1.csv"))

df = data.frame(biasMat1)
colnames(df) <- c('tas' ,'pr', 'rsds','rlds' , 'hurs')
write.table(df, paste0(wd, "/bias_stats1.csv"))
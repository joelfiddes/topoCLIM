require(qmap)
require(ncdf4)
require(caTools)
require(PCICt) # converts to 365 day clim model calender
require(pracma)
require(hydroGOF)
require(wesanderson)
require(stargazer)
# Setup ========================================================================

# POI (WFJ) - to extract cordex gridbox
mylon = 9.80490
mylat = 46.83066

#units
#		pr:standard_name = "precipitation_flux" ;
#		pr:long_name = "Precipitation" ;
#		pr:units = "kg m-2 s-1" 

# Cordex data
hist_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/historical/r1i1p1/REMO2015/v1/day/pr/pr_EUR-22_MOHC-HadGEM2-ES_historical_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"
rcp26_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/rcp26/r1i1p1/REMO2015/v1/day/pr/pr_EUR-22_MOHC-HadGEM2-ES_rcp26_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"
rcp85_file= "/home/joel/sim/qmap/CORDEX/output/EUR-22/GERICS/MOHC-HadGEM2-ES/rcp85/r1i1p1/REMO2015/v1/day/pr/pr_EUR-22_MOHC-HadGEM2-ES_rcp85_r1i1p1_GERICS-REMO2015_v1_day_TS.nc"

# Obs Era5 data downscaled by toposcale, can be resampled

obsfile = "/home/joel/sim/qmap/wfj_long_no_plapse.csv"

# Obs 1980-01-01 - 2017-12-31 ==========================================================================
#read obs variable TA
obs=read.csv(obsfile)
pr_obs=obs$PINT				# !!HARDCODE!!

# aggregate obs to daily resolution
days=substr(obs$datetime,1,10)
d =aggregate(pr_obs, list(days), 'mean')
obs_daily=d$x
obs_daily_date=d$Group.1



# Cordex historical "1970-01-01 12" UTC" to 2005-12-30 12"======================
modfile=hist_file
# define variable and gridbox
nc=nc_open(modfile)
pr =ncvar_get(nc, 'pr') 		# !!HARDCODE!!
lon =ncvar_get(nc, 'lon')
lat =ncvar_get(nc, 'lat')
myx =which.min(abs(lon - mylon)) 
myy = which.min(abs(lat - mylat)) 

# extract timeseries corresponding to qmap (all)
pr_all = pr[myx,myy,]*60*60


# define time
time = ncvar_get( nc,'time')
time2=time
z <- time2*24*60*60 #convert days to seconds
#origin=substr(nc$dim$time$units,13,20)
origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section


cal <- "360_day"
seconds.per.day <- 86400
ts.dat.days <- time
origin.pcict <- as.PCICt(origin, cal)
ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
datesPl <- ts.dat.pcict 

# convert 360 to 365/366 calender with interpolation
#greg_cal = seq(as.Date('1970-01-01'), as.Date('2005-12-31'), "days") # hardcoded
greg_cal_hist = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10))+1, "days") # add plus 1 here as cordex always ends on dec 30

# get index for interp1d
greg_unix = as.numeric(as.POSIXct(greg_cal_hist, format="%Y-%m-%d"))
cordCal_unix = as.numeric(as.POSIXct( as.Date(substr(datesPl,1,10)), format="%Y-%m-%d"))

# NAs generated in cordCal_unix as dates are made in the pcict calender conversion that dont exists ie Feb 30
# we remove these values completely representing 0.48% of values (63/12960)
vals2rm = which(is.na(cordCal_unix)==T)
vals2keep = which(!is.na(cordCal_unix)==T)


overlap_index = which(greg_unix %in% cordCal_unix)
pr_365_minus1 = interp1(x=cordCal_unix[vals2keep], y=pr_all[vals2keep], xi = greg_unix[1:(length(greg_unix)-1)], method = c("linear"))

# we miss very last dec 31 as x only goes to dec 30 and inter1 cannot interp beyond range of x (extrapolate). thats why we do this -1: xi = greg_unix[1:(length(greg_unix)-1)]
# to fix we simply repeat last value
pr_365 = c(pr_365_minus1, pr_365_minus1[length(pr_365_minus1)])

#datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
days_cordex=greg_cal_hist
years_cordex=as.numeric(substr(greg_cal_hist,1,4))
start_cp_cord = which(days_cordex==days[1])
end_cp_cord=which(days_cordex==days[length(days)])

if(length(end_cp_cord)==0){ # if true this mean that obs extend beyond end_cp_cord of the cordex data - this is the case for cordex historical
end_cp_cord= length(days_cordex) # then simply take last cordex date
}
if(length(start_cp_cord)==0){ # if true this mean that obs starts before cordex data - this is the case for cordex rcp data
start_cp_cord = 1 # then simply take first cordex date
}
# extract timeseries corresponding to obs (subset in time/space)
hist_cp=pr_365[start_cp_cord:end_cp_cord]

# define common period
start_cp= greg_cal_hist[start_cp_cord]
end_cp = greg_cal_hist[end_cp_cord]

# cut obs to common period cp
start_obs = which(obs_daily_date==start_cp)
end_obs=which(obs_daily_date==end_cp)
obs_daily_cp = obs_daily[start_obs: end_obs]

# get gmap pars
pars = fitQmap(obs_daily_cp,hist_cp, method = "QUANT")

# seasonal bias correction
month_cp = substring(seq(start_cp, end_cp, "days"),6,7)
summer_index=which(month_cp=="04"| month_cp=="05"| month_cp=="06"| month_cp=="07"| month_cp=="08"| month_cp=="09")
winter_index=which(month_cp=="10"| month_cp=="11"| month_cp=="12"| month_cp=="01"| month_cp=="02"| month_cp=="03")
pars_summer = fitQmap(obs_daily_cp[summer_index],hist_cp[summer_index], method = "QUANT")
pars_winter = fitQmap(obs_daily_cp[winter_index],hist_cp[winter_index], method = "QUANT")



# do qmap
hist = doQmap(pr_365,pars)

# do seasonal qmap
month_cp = substring(greg_cal_hist,6,7)
summer_index=which(month_cp=="04"| month_cp=="05"| month_cp=="06"| month_cp=="07"| month_cp=="08"| month_cp=="09")
winter_index=which(month_cp=="10"| month_cp=="11"| month_cp=="12"| month_cp=="01"| month_cp=="02"| month_cp=="03")
hist_summer = doQmap(pr_365[summer_index],pars_summer)
hist_winter = doQmap(pr_365[winter_index],pars_winter)
hist_season<-1:length(hist)
hist_season[summer_index]<-hist_summer
hist_season[winter_index]<-hist_winter
#annual aggregated timeseries
hist_year=aggregate(hist, list(years_cordex), 'mean')
hist_season_year=aggregate(hist_season, list(years_cordex), 'mean')
pr_365_year=aggregate(pr_365, list(years_cordex), 'mean')
obs_year=aggregate(obs_daily, list(substring(obs_daily_date,1,4)), 'mean')

# rmove invalid first and last datapoint (incomplete years)
#hist_year$x[1]<-NA
#hist_year$x[length(hist_year$x)]<-NA
plot(hist_year, type='l')




# Cordex rcp26 =======================================================================
modfile=rcp26_file
# define variable and gridbox
nc=nc_open(modfile)
pr =ncvar_get(nc, 'pr') 		# !!HARDCODE!!
lon =ncvar_get(nc, 'lon')
lat =ncvar_get(nc, 'lat')
myx =which.min(abs(lon - mylon)) 
myy = which.min(abs(lat - mylat)) 

# extract timeseries corresponding to qmap (all)
pr_all = pr[myx,myy,]*60*60 # mm/s >mm/h


# define time
time = ncvar_get( nc,'time')
time2=time
z <- time2*24*60*60 #convert days to seconds
#origin=substr(nc$dim$time$units,13,20)
origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section


cal <- "360_day"
seconds.per.day <- 86400
ts.dat.days <- time
origin.pcict <- as.PCICt(origin, cal)
ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
datesPl <- ts.dat.pcict 

# convert 360 to 365/366 calender with interpolation
greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10))+1, "days") # add plus 1 here as cordex always ends on dec 30

# get index for interp1d
greg_unix = as.numeric(as.POSIXct(greg_cal_rcp, format="%Y-%m-%d"))
cordCal_unix = as.numeric(as.POSIXct( as.Date(substr(datesPl,1,10)), format="%Y-%m-%d"))

# NAs generated in cordCal_unix as dates are made in the pcict calender conversion that dont exists ie Feb 30
# we remove these values completely representing 0.48% of values (63/12960)
vals2rm = which(is.na(cordCal_unix)==T)
vals2keep = which(!is.na(cordCal_unix)==T)


overlap_index = which(greg_unix %in% cordCal_unix)
pr_365_minus1 = interp1(x=cordCal_unix[vals2keep], y=pr_all[vals2keep], xi = greg_unix[1:(length(greg_unix)-1)], method = c("linear"))

# we miss very last dec 31 as x only goes to dec 30 and inter1 cannot interp beyond range of x (extrapolate). thats why we do this -1: xi = greg_unix[1:(length(greg_unix)-1)]
# to fix we simply repeat last value
pr_365 = c(pr_365_minus1, pr_365_minus1[length(pr_365_minus1)])

#datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
days_cordex=greg_cal_rcp
years_cordex=as.numeric(substr(greg_cal_rcp,1,4))
#start_cp_cord = which(days_cordex==days[1])
#end_cp_cord=which(days_cordex==days[length(days)])

#if(length(end_cp_cord)==0){ # if true this mean that obs extend beyond end_cp_cord of the cordex data - this is the case for cordex historical
#end_cp_cord = length(days_cordex) # then simply take last cordex date
#}
#if(length(start_cp_cord)==0){ # if true this mean that obs starts before cordex data - this is the case for cordex rcp data
#start_cp_cord = 1 # then simply take first cordex date
#}
## extract timeseries corresponding to obs (subset in time/space)
#hist_cp=pr_365[start_cp_cord:end_cp_cord]

## define common period
#start_cp= greg_cal_rcp[start_cp_cord]
#end_cp = greg_cal_rcp[end_cp_cord]

## cut obs to common period cp
#start_obs = which(obs_daily_date==start_cp)
#end_obs=which(obs_daily_date==end_cp)
#obs_daily_cp = obs_daily[start_obs: end_obs]

# get gmap pars
#pars = fitQmap(obs_daily_cp,hist_cp, method = "QUANT")

# do qmap
rcp26 = doQmap(pr_365,pars)

#annual aggregated timeseries
rcp26_year=aggregate(rcp26, list(years_cordex), 'mean')

# rmove invalid first and last datapoint (incomplete years)
#hist_year$x[1]<-NA
#hist_year$x[length(hist_year$x)]<-NA





# Cordex rcp85 2006-2100=======================================================================
modfile=rcp85_file
# define variable and gridbox
nc=nc_open(modfile)
pr =ncvar_get(nc, 'pr') 		# !!HARDCODE!!
lon =ncvar_get(nc, 'lon')
lat =ncvar_get(nc, 'lat')
myx =which.min(abs(lon - mylon)) 
myy = which.min(abs(lat - mylat)) 

# extract timeseries corresponding to qmap (all)
pr_all = pr[myx,myy,]*60*60


# define time
time = ncvar_get( nc,'time')
time2=time
z <- time2*24*60*60 #convert days to seconds
#origin=substr(nc$dim$time$units,13,20)
origin = unlist(strsplit(nc$dim$time$units, " "))[[3]]
origin2 = substring(origin, 1,10) # extract just yyyy-mm-dd section


cal <- "360_day"
seconds.per.day <- 86400
ts.dat.days <- time
origin.pcict <- as.PCICt(origin, cal)
ts.dat.pcict <- origin.pcict + (ts.dat.days * seconds.per.day)
datesPl <- ts.dat.pcict 

# convert 360 to 365/366 calender with interpolation
greg_cal_rcp = seq(as.Date(substring(datesPl[1],1,10)), as.Date(substring(datesPl[length(datesPl)],1,10))+1, "days") # add plus 1 here as cordex always ends on dec 30

# get index for interp1d
greg_unix = as.numeric(as.POSIXct(greg_cal_rcp, format="%Y-%m-%d"))
cordCal_unix = as.numeric(as.POSIXct( as.Date(substr(datesPl,1,10)), format="%Y-%m-%d"))

# NAs generated in cordCal_unix as dates are made in the pcict calender conversion that dont exists ie Feb 30
# we remove these values completely representing 0.48% of values (63/12960)
vals2rm = which(is.na(cordCal_unix)==T)
vals2keep = which(!is.na(cordCal_unix)==T)


overlap_index = which(greg_unix %in% cordCal_unix)
pr_365_minus1 = interp1(x=cordCal_unix[vals2keep], y=pr_all[vals2keep], xi = greg_unix[1:(length(greg_unix)-1)], method = c("linear"))

# we miss very last dec 31 as x only goes to dec 30 and inter1 cannot interp beyond range of x (extrapolate). thats why we do this -1: xi = greg_unix[1:(length(greg_unix)-1)]
# to fix we simply repeat last value
pr_365 = c(pr_365_minus1, pr_365_minus1[length(pr_365_minus1)])

#datesPl<-ISOdatetime(origin2,0,0,0,0,0,tz='UTC') + z #dates sequence
days_cordex=greg_cal_rcp
years_cordex=as.numeric(substr(greg_cal_rcp,1,4))
#start_cp_cord = which(days_cordex==days[1])
#end_cp_cord=which(days_cordex==days[length(days)])

#if(length(end_cp_cord)==0){ # if true this mean that obs extend beyond end_cp_cord of the cordex data - this is the case for cordex historical
#end_cp_cord = length(days_cordex) # then simply take last cordex date
#}
#if(length(start_cp_cord)==0){ # if true this mean that obs starts before cordex data - this is the case for cordex rcp data
#start_cp_cord = 1 # then simply take first cordex date
#}
## extract timeseries corresponding to obs (subset in time/space)
#hist_cp=pr_365[start_cp_cord:end_cp_cord]

## define common period
#start_cp= greg_cal_rcp[start_cp_cord]
#end_cp = greg_cal_rcp[end_cp_cord]

## cut obs to common period cp
#start_obs = which(obs_daily_date==start_cp)
#end_obs=which(obs_daily_date==end_cp)
#obs_daily_cp = obs_daily[start_obs: end_obs]

# get gmap pars
#pars = fitQmap(obs_daily_cp,hist_cp, method = "QUANT")

# do qmap
rcp85 = doQmap(pr_365,pars)

#annual aggregated timeseries
rcp85_year=aggregate(rcp85, list(years_cordex), 'mean')

# rmove invalid first and last datapoint (incomplete years)
#hist_year$x[1]<-NA
#hist_year$x[length(hist_year$x)]<-NA
plot(rcp85_year, type='l')

# summary plot=======================================================================
mycol = wes_palette("Zissou1", 3, type = "continuous")

lwd=3
plot(rcp26_year$Group.1,rcp26_year$x*24*365, xlim=c(1970,2100),ylim=c(1000,3000), type='l', col=mycol[1], lwd=lwd, main="Mean annual air temp Weissfluhjoch", ylab="Air temperacture (K)", xlab=" ")
lines(rcp85_year$Group.1,rcp85_year$x*24*365,  type='l', col=mycol[2],lwd=lwd)
lines(hist_year$Group.1,hist_year$x*24*365, , type='l', col=mycol[3],lwd=lwd)
lines(pr_365_year$Group.1,pr_365_year$x*24*365, type='l', col='grey',lwd=lwd,lty=2)
lines(obs_year$Group.1, obs_year$x*24*365, lty=2, lwd=lwd)
legend("bottomright", c("Hist", "RCP2.6", "RCP8.5", "Hist_noqmap", "obs"), col=c(mycol[3],mycol[1],mycol[2], mycol[3], 'black'),lwd=3, lty=c(1,1,1,1,2))


#plot(datesPl,rcp26_all, type='l', col='blue', lwd=lwd, main="Mean annual air temp Weissfluhjoch", ylab="Air temperacture (K)", xlab=" ")
#lines(rcp85_all, type='l', col='red',lwd=lwd)
#lines(pr_all, ylim=c(270,280), type='l', col='green',lwd=lwd)
#legend("bottomright", c("Hist", "RCP2.6", "RCP8.5"), col=c("green","blue", "red"),lwd=3)







# Evaluation historical (1980-01-01 - 2005-12-31)
#===============================================================================
# prepare data
#===============================================================================
obs = obs_daily[1:length(hist_cp)]
hist_noqmap = hist_cp
hist_qmap =hist[start_cp_cord:end_cp_cord]
hist_season_qmap =hist_season[start_cp_cord:end_cp_cord]
dates_cp = greg_cal_hist[start_cp_cord:end_cp_cord]
dates_month = substring(dates_cp,1,7)
obs_month = aggregate(obs, list(dates_month), mean)[,2] 
hist_month = aggregate(hist_noqmap, list(dates_month), mean)[,2] 
hist_qmap_month = aggregate(hist_qmap, list(dates_month), mean)[,2]
hist_season_qmap_month = aggregate(hist_season_qmap, list(dates_month), mean)[,2]

#===============================================================================
# summary stats
#===============================================================================
cor1 =cor(obs_month, hist_month)
cor2 =cor(obs_month, hist_qmap_month)
cor3 =cor(obs_month, hist_season_qmap_month)

rms1 =rmse(obs_month, hist_month)
rms2 =rmse(obs_month, hist_qmap_month)
rms3 =rmse(obs_month, hist_season_qmap_month)

m1=mean(obs_month)
m2=mean(hist_month)
m3=mean(hist_qmap_month)
m3=mean(hist_season_qmap_month)
sumtab=round(as.data.frame(rbind( c(m1,m2,m3,m3), c(NA, cor1,cor2,cor3), c(NA,rms1,rms2,rms3))),2)
row.names(sumtab) <- c("MEAN", "R-COR" ,"RMSE")
names(sumtab) <-c("OBS","Hist", "Hist_QMAP", "Hist_QMAPSEASON")
stargazer(sumtab, summary=FALSE)

#% Table created by stargazer v.5.2.2 by Marek Hlavac, Harvard University. E-mail: hlavac at fas.harvard.edu
#% Date and time: Fr, Jul 03, 2020 - 12:27:21
#\begin{table}[!htbp] \centering 
#  \caption{} 
#  \label{} 
#\begin{tabular}{@{\extracolsep{5pt}} ccccc} 
#\\[-1.8ex]\hline 
#\hline \\[-1.8ex] 
# & OBS & Hist & Hist\_QMAP & Hist\_QMAPSEASON \\ 
#\hline \\[-1.8ex] 
#MEAN & $0.170$ & $0.240$ & $0.170$ & $0.170$ \\ 
#R-COR & $$ & $$-$0.090$ & $$-$0.100$ & $0.210$ \\ 
#RMSE & $$ & $0.180$ & $0.130$ & $0.110$ \\ 
#\hline \\[-1.8ex] 
#\end{tabular} 
#\end{table}

#===============================================================================
# scatter plot
#===============================================================================
pdf("/home/joel/manuscripts/qmap/plots/P_SCATTER.pdf")
lwd=3
plot(obs_month, hist_month,col=mycol[2], ylab="Modelled Precipitaion (kg m-2 s-1)", xlab="Observed Precipitation (kg m-2 s-1)", ylim=c(0,0.7), xlim=c(0,0.7))
points(obs_month, hist_qmap_month,col=mycol[3])
points(obs_month, hist_season_qmap_month,col=mycol[1])
abline(0,1)
abline(lm(hist_month~obs_month),col=mycol[2], lwd=lwd)
abline(lm(hist_qmap_month~obs_month),col=mycol[3], lwd=lwd)
abline(lm(hist_season_qmap_month~obs_month),col=mycol[1], lwd=lwd)
legend("topright", c("1:1 line", "CLIM", "CLIM_QM","CLIM_QMSEASON"), col=c("black",mycol[2], mycol[3],  col=mycol[1]),lty=c(2,1,1,1) ,lwd=lwd)
dev.off()




#===============================================================================
# ecdf
#===============================================================================
pdf("/home/joel/manuscripts/qmap/plots/P_CDF.pdf")
lwd=2
 plot(ecdf(obs_month), ylab="cdf", xlab="P (kg m-2 s-1)", main=" " , lwd=lwd,lty=1)
 lines(ecdf(hist_month),  col=mycol[2], lwd=lwd)
 lines(ecdf(hist_qmap_month),  col=mycol[3], lwd=lwd)
 lines(ecdf(hist_season_qmap_month),  col=mycol[1], lwd=lwd)
legend("bottomright", c("OBS", "CLIM", "CLIM_QM","CLIM_QMSEASON"), col=c("black",mycol[2], mycol[3],  col=mycol[1]),lty=c(1,1,1,1) ,lwd=lwd)
dev.off()

# daily ecdf
 plot(ecdf(obs), ylab="cdf", xlab="TA (K)", main=" ", lwd=lwd)
 lines(ecdf(hist_noqmap), col='red', lwd=lwd)
 lines(ecdf(hist_qmap), col='blue', lwd=lwd)
legend("bottomright", c("OBS", "CLIM", "CLIM_QM"), col=c("black","red", "blue"), lwd=lwd)
# day plot
doy=substr(unique(dates_cp), 6,10)


#===============================================================================
# DOY
#===============================================================================
pdf("/home/joel/manuscripts/qmap/plots/P_DOY.pdf")
doy_obs = aggregate(obs, list(doy), mean)
doy_clim = aggregate(hist_noqmap, list(doy), mean)
doy_qm = aggregate(hist_qmap, list(doy), mean)
doy_season_qm = aggregate(hist_season_qmap, list(doy), mean)
#
plot(runmean(doy_obs$x,30), type='l', ylim=c(0,0.4),lwd=lwd,lty=2, xlab='DOY', ylab='P (kg m-2 s-1)')
lines(runmean(doy_clim$x,30), type='l', col=mycol[2], lwd=lwd)
lines(runmean(doy_qm$x,30), type='l', col=mycol[3],lwd=lwd)
lines(runmean(doy_season_qm$x,30),col=mycol[1], type='l',lwd=3)
legend("topright", c("OBS", "CLIM", "CLIM_QM","CLIM_QMSEASON"), col=c("black",mycol[2], mycol[3],  col=mycol[1]),lty=c(2,1,1,1) ,lwd=lwd)
dev.off()




# we evaluate with rcp data in period 2006-2017























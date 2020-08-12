# Disaggregate daily to sub daily timeseries in order to drive models


windDir <- function(u, v) {
  (180 / pi) * atan(u/v) + ifelse(v>0,180,ifelse(u>0,360,0))
}
#source: http://stackoverflow.com/questions/8673137/calculating-wind-direction-from-u-and-v-components-of-the-wind-using-lapply-or-i

################################################################################################################
#					calc wind speed from vectors
################################################################################################################
#function calc wind speed [liston elder 2006]
windSpd<-function(u,v){

ws<-sqrt(u^2+v^2)
#ws<-180/pi*s
return(ws)
}



obsfile = "/home/joel/sim/qmap/wfj_long.csv"#
obs=read.csv(obsfile)

 days=substr(obs$datetime,1,10)
 obs_day =aggregate(obs, list(days), 'mean')

obs_day_dat <- obs_day[,c(12,9,7,6,4,5,13,3)]

# conversions
obs_day_dat$PINT<- obs_day_dat$PINT/3600
obs_day_dat$RH<- obs_day_dat$RH*100

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_tas')
hist_tas <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_hurs')
hist_hurs <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_pr')
hist_pr <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_ps')
hist_ps <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_rlds')
hist_rlds <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_rsds')
hist_rsds <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_uas')
hist_uas <- df

load('/home/joel/sim/qmap/test/pyout/aresult/aqmap_results/hist_vas')
hist_vas <- df

allModels <- unique(c(names(hist_tas),
names(hist_pr),
names(hist_ps),
names(hist_rlds),
names(hist_rsds),
names(hist_uas),
names(hist_vas)))

allModels <- allModels[- which(allModels=='Date')]


for (model in allModels){
	print(model)
idtas = which(names(hist_tas)==model)
idhurs = which(names(hist_hurs)==model)
idpr = which(names(hist_pr)==model)
idps = which(names(hist_ps)==model)
idrlds = which(names(hist_rlds)==model)
idrsds = which(names(hist_rsds)==model)
iduas = which(names(hist_uas)==model)
idvas = which(names(hist_vas)==model)


tas <-hist_tas[,idtas]
hurs <-hist_hurs[,idhurs]
pr <-hist_pr[,idpr]
ps <-hist_ps[,idps]
rlds <-hist_rlds[,idrlds]
rsds <-hist_rsds[,idrsds]
uas <-hist_uas[,iduas]
vas <-hist_vas[,idvas]

ws <- windSpd(uas,vas)
wd <- windDir(uas,vas)

df=data.frame(tas,hurs, pr,ps,rlds,rsds,ws,wd)

if( dim(df)[2]!=8){print(paste0("Not enough vars ", model)); next}
dfscale=scale(df)
normFact = colMeans(as.matrix(obs_day_dat))
for (i in dim(df)[1]){
	a=  which.min(abs(rowSums(sweep((as.matrix(obs_day_dat)/ colMeans(as.matrix(obs_day_dat))),2,(as.numeric(df[i,])/ colMeans(as.matrix(obs_day_dat)))), na.rm=T)))
	obs_day_dat[a,]
	df[i,]
	}



#assign( paste0("DAT_",gsub('-','',model)), df )
}
# construct obs timeseries


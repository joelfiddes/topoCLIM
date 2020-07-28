# decompose daily to subdaily using analog approach
# caveats subdaily distribution may change yadahyadah
# consistent as daily only freq with all 8 variable
# more efficient download etc easy to handle

require(ncdf4)
# POI (WFJ) - to extract cordex gridbox
mylon = 9.80490
mylat = 46.83066




vars=c('tas', 'pr', 'ps', 'hurs', 'rsds','rlds', 'uas', 'vas')

# constrct test climate data CNRM (most complete)
nc = nc_open("/home/joel/sim/qmap/CORDEX/output/EUR-44/final_results/_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_CNRM-ALADIN53_v1_day__TS.nc_5VARS.nc")
nc = nc_open("/home/joel/sim/qmap/CORDEX/output/EUR-44/final_results/_CNRM-CERFACS-CNRM-CM5_historical_r1i1p1_HMS-A_TS.nc_6VARS.nc")


lon =ncvar_get(nc, 'lon')
lat =ncvar_get(nc, 'lat')
myx =which.min(abs(lon - mylon))
myy = which.min(abs(lat - mylat))
tas=ncvar_get(nc, 'tas')
hurs=ncvar_get(nc, 'hurs')
ps=ncvar_get(nc, 'ps')
pr=ncvar_get(nc, 'pr')
rlds=ncvar_get(nc, 'rlds')
rsds=ncvar_get(nc, 'rsds')

# construct dataset for cordex
cordex_df =data.frame(rlds[myx, myy, ], rsds[myx, myy, ], ps[myx, myy, ]*100, pr[myx, myy, ]*3600, hurs[myx, myy, ]/100, tas[myx, myy, ])





# subdaily data (ERA5)
obsfile = "/home/joel/sim/qmap/wfj_long.csv"
obs=read.csv(obsfile)
> names(obs)

 [1] "datetime" "DW"       "ILWR"     "ISWR"     "P"        "PINT"
 [7] "PSUM"     "RH"       "Rf"       "Sf"       "TA"       "VW"

 # aggregate obs to daily resolution
 days=substr(obs$datetime,1,10)
 d =aggregate(obs, list(days), 'mean')
 #obs_daily=d$x
 #obs_daily_date=d$Group.1

 # comparable
 obs_df = d[c(4,5,6,7,9,12)]


for (i in dim(cordex_df)[1]){

a=  rowSums(sweep(as.matrix(obs_df),2,unname(as.numeric(cordex_df[i,]))))

}
b[,3:5] <- t(apply(b[,3:5], 1, function(x) x-c))

args = commandArgs(trailingOnly=TRUE)
root = args[1]
sample = args[2]
indir=paste0(root, '/s',sample)



fsmdir=paste0(indir,"/fsm/")
dir.create(fsmdir)

results1 = list.files(path = paste0(indir,"/aqmap_results/"),full.names=T )
results = results1[ grep('_season_',results1)]

print(results)
hist = results[grep("hist",results)]
rcp26 = results[grep("rcp26",results)]
rcp85 = results[grep("rcp85",results)]

#DW

typ='HIST'
load(hist[1] )
RH<-df
load(hist[2] )
PINT<-df
load(hist[3] )
P<-df
load(hist[4] )
ILWR<-df
load(hist[5] )
ISWR<-df
load(hist[6] )
TA<-df
load(hist[7] )
TAMAX<-df
load(hist[8] )
TAMIN<-df
load(hist[9] )
VW<-df[1:length(names(df))-1]



#DW = (180 / pi) * atan(uas/vas) + ifelse(vas>0,180,ifelse(uas>0,360,0))
#VW = sqrt(uas^2+vas^2)

models = names(df)[1:length(names(df))-1]
dates = df[length(df)]
for (model in models){
	print(model)
	modellist=list()
	modellist['datetime'] <- dates
	#modellist['DW']<- DW[model]
	modellist['ILWR']<- ILWR[model]
	modellist['ISWR']<- ISWR[model] # qmap can make non-zero night time rad which is unphysical
	modellist['P']<- P[model]
	modellist['PINT']<- PINT[model]
	modellist['RH']<- RH[model]
	modellist['TA']<- TA[model]
	modellist['VW']<- VW[model]
	modellist['TAMIN']<- TAMIN[model]
	modellist['TAMAX']<- TAMAX[model]
	df = as.data.frame(modellist)
	df[c(2:4, 6:10)] =round(df[c(2:4, 6:10)],1)
	df[5 ] = round(df[5],6)
	write.csv(df, paste0(fsmdir, model, "_",typ,"_Q.txt"), row.names=FALSE)
	
}


typ='RCP26'
load(rcp26[1] )
RH<-df
load(rcp26[2] )
PINT<-df
load(rcp26[3] )
P<-df
load(rcp26[4] )
ILWR<-df
load(rcp26[5] )
ISWR<-df
load(rcp26[6] )
TA<-df
load(rcp26[7] )
TAMAX<-df
load(rcp26[8] )
TAMIN<-df
load(rcp26[9] )
VW<-df[1:length(names(df))-1]



#DW = (180 / pi) * atan(uas/vas) + ifelse(vas>0,180,ifelse(uas>0,360,0))
#VW = sqrt(uas^2+vas^2)

models = names(df)[1:length(names(df))-1]
dates = df[length(df)]
for (model in models){
	print(model)
	modellist=list()
	modellist['datetime'] <- dates
	#modellist['DW']<- DW[model]
	modellist['ILWR']<- ILWR[model]
	modellist['ISWR']<- ISWR[model] # qmap can make non-zero night time rad which is unphysical
	modellist['P']<- P[model]
	modellist['PINT']<- PINT[model]
	modellist['RH']<- RH[model]
	modellist['TA']<- TA[model]
	modellist['VW']<- VW[model]
	modellist['TAMIN']<- TAMIN[model]
	modellist['TAMAX']<- TAMAX[model]
	df = as.data.frame(modellist)
	df[c(2:4, 6:10)] =round(df[c(2:4, 6:10)],1)
	df[5 ] = round(df[5],6)
	write.csv(df, paste0(fsmdir, model, "_",typ,"_Q.txt"), row.names=FALSE)
	
}

typ='RCP85'
load(rcp85[1] )
RH<-df
load(rcp85[2] )
PINT<-df
load(rcp85[3] )
P<-df
load(rcp85[4] )
ILWR<-df
load(rcp85[5] )
ISWR<-df
load(rcp85[6] )
TA<-df
load(rcp85[7] )
TAMAX<-df
load(rcp85[8] )
TAMIN<-df
load(rcp85[9] )
VW<-df[1:length(names(df))-1]



#DW = (180 / pi) * atan(uas/vas) + ifelse(vas>0,180,ifelse(uas>0,360,0))
#VW = sqrt(uas^2+vas^2)

models = names(df)[1:length(names(df))-1]
dates = df[length(df)]
for (model in models){
	print(model)
	modellist=list()
	modellist['datetime'] <- dates
	#modellist['DW']<- DW[model]
	modellist['ILWR']<- ILWR[model]
	modellist['ISWR']<- ISWR[model] # qmap can make non-zero night time rad which is unphysical
	modellist['P']<- P[model]
	modellist['PINT']<- PINT[model]
	modellist['RH']<- RH[model]
	modellist['TA']<- TA[model]
	modellist['VW']<- VW[model]
	modellist['TAMIN']<- TAMIN[model]
	modellist['TAMAX']<- TAMAX[model]
	df = as.data.frame(modellist)
	df[c(2:4, 6:10)] =round(df[c(2:4, 6:10)],1)
	df[5 ] = round(df[5],6)
	write.csv(df, paste0(fsmdir, model, "_",typ,"_Q.txt"), row.names=FALSE)
	
}
require(raster)
require(wesanderson)
#https://github.com/karthik/wesanderson/blob/master/R/colors.R
Zissou1 = c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
pal <- Zissou1
discrete = pal[1:3]

wesblueT=c("#3B9AB280") 
wesgreenT=c("#EBCC2A80") 
wesredT=c("#F21A0080") 

wesblue=c("#3B9AB2") 
wesgreen=c("#EBCC2A") 
wesred=c("#F21A00") 

mycol = wes_palette("Zissou1", 3, type = "continuous")
glacier =shapefile("/home/joel/data/glims/alps/rgi/11_rgi60_CentralEurope.shp")
ele = raster("/home/joel/sim/qmap/ch_tmapp2/eleCH.tif")
sweHIST = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/HIST_2__map.tif")


swercp26near = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP26_2_2030-01-01_2050-12-31__map.tif")
swercp26far = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP26_2_2080-01-01_2099-12-31__map.tif")


swercp85near = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP85_2_2030-01-01_2050-12-31__map.tif")
swercp85far = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP85_2_2080-01-01_2099-12-31__map.tif")


root ="/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP"





#=====================================================================================
# NEW CODE
#=====================================================================================

ele = raster("/home/joel/sim/qmap/ch_tmapp2/eleCH.tif")


slim=1
#=====================================================================================
# Hist
#=====================================================================================
step=50
llvec=seq(400,4800,step)
slim=4
el=llvec[1:length(llvec)-1]

swe_raster = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/HIST_2__map.tif")
swe_raster[swe_raster>slim]<-NA

meanlist=list()
sdlist= list()

for (i in 1:(length(llvec)-1)){

ll1=100
ll=llvec[i]
ul=llvec[i+1]

m <- c(ll1, ll,NA,ll, ul, 1,  ul, llvec[length(llvec)], NA  )
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(ele, rclmat, include.lowest=TRUE)
b9 <- rc*swe_raster

meanlist[[i]] = cellStats(b9, mean,na.rm=T)
sdlist[[i]]<- cellStats(b9,sd, na.rm=T)

}

ml1 = unlist(meanlist)
sl1 = unlist(sdlist)



#=====================================================================================
# rcp85 near
#=====================================================================================
swe_raster = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP85_2_2030-01-01_2050-12-31__map.tif")
swe_raster[swe_raster>slim]<-NA

meanlist=list()
sdlist= list()

for (i in 1:(length(llvec)-1)){

ll1=100
ll=llvec[i]
ul=llvec[i+1]


m <- c(ll1, ll,NA,ll, ul, 1,  ul, llvec[length(llvec)], NA  )
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(ele, rclmat, include.lowest=TRUE)
b9 <- rc*swe_raster

meanlist[[i]] = cellStats(b9, mean,na.rm=T)
sdlist[[i]]<- cellStats(b9,sd,na.rm=T)

}

ml85nr = unlist(meanlist)
sl85nr = unlist(sdlist)

#=====================================================================================
# rcp85 far
#=====================================================================================
swe_raster = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP85_2_2080-01-01_2099-12-31__map.tif")
swe_raster[swe_raster>slim]<-NA

meanlist=list()
sdlist= list()

for (i in 1:(length(llvec)-1)){

ll1=100
ll=llvec[i]
ul=llvec[i+1]


m <- c(ll1, ll,NA,ll, ul, 1,  ul, llvec[length(llvec)], NA  )
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(ele, rclmat, include.lowest=TRUE)
b9 <- rc*swe_raster

meanlist[[i]] = cellStats(b9, mean,na.rm=T)
sdlist[[i]]<- cellStats(b9,sd,na.rm=T)

}

ml85fr = unlist(meanlist)
sl85fr = unlist(sdlist)


#=====================================================================================
# rcp26 nr
#=====================================================================================
swe_raster = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP26_2_2030-01-01_2050-12-31__map.tif")
swe_raster[swe_raster>slim]<-NA

meanlist=list()
sdlist= list()

for (i in 1:(length(llvec)-1)){

ll1=100
ll=llvec[i]
ul=llvec[i+1]


m <- c(ll1, ll,NA,ll, ul, 1,  ul, llvec[length(llvec)], NA  )
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(ele, rclmat, include.lowest=TRUE)
b9 <- rc*swe_raster

meanlist[[i]] = cellStats(b9, mean,na.rm=T)
sdlist[[i]]<- cellStats(b9,sd,na.rm=T)

}

ml26nr = unlist(meanlist)
sl26nr = unlist(sdlist)

#=====================================================================================
# rcp26 far
#=====================================================================================
swe_raster = raster("/home/joel/sim/qmap/PAPER_RESULTS/ch_tmapp_50-HYP/RCP26_2_2080-01-01_2099-12-31__map.tif")
swe_raster[swe_raster>slim]<-NA

meanlist=list()
sdlist= list()

for (i in 1:(length(llvec)-1)){

ll1=100
ll=llvec[i]
ul=llvec[i+1]


m <- c(ll1, ll,NA,ll, ul, 1,  ul, llvec[length(llvec)], NA  )
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(ele, rclmat, include.lowest=TRUE)
b9 <- rc*swe_raster

meanlist[[i]] = cellStats(b9, mean,na.rm=T)
sdlist[[i]]<- cellStats(b9,sd,na.rm=T)

}

ml26fr = unlist(meanlist)
sl26fr = unlist(sdlist)

pdf("/home/joel/manuscripts/qmap/plots/hypsometry_HS.pdf")
par(mfrow=c(1,2))
lwd=2


#rcp26
mymean = ml1
#ml1[74:88]<-1.186
#ml1[is.na(ml1)] <- rev(na.omit(ml1))[1] 

sl1[is.na(sl1)]<-0 
mysd = sl1
plot(mymean,el,lwd=lwd, xlim=c(0,3), ylim=c(2000,5000), typ="l", col=wesgreen, xlab="Mean snow depth (m)", ylab="Elevation (m asl)", main="RCP26")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=wesgreenT,border='NA')

mymean = ml26nr
mysd = sl26nr
lines(mymean,el,lwd=lwd,col=wesblue)
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=wesblueT,border='NA')


mymean = ml26fr
mysd = sl26fr
lines(mymean,el,lwd=lwd,col=wesred)
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=wesredT,border='NA')
legend("bottomright", legend=c("1980-2010", "2030-60", "2070-99"), col=c(wesgreen,wesblue,wesred),lty=1)
# rcp 85
mymean = ml1
mysd = sl1
plot(mymean,el,lwd=lwd,  xlim=c(0,3), ylim=c(2000,5000), typ="l", col=wesgreen, xlab="Mean snow depth (m)", ylab="Elevation (m asl)", main="RCP85")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=wesgreenT,border='NA')

mymean = ml85nr
mysd = sl85nr
lines(mymean,el,lwd=lwd, col=wesblue)
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=wesblueT,border='NA')


mymean = ml85fr
mysd = sl85fr
lines(mymean,el,lwd=lwd, col=wesred)
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=wesredT,border='NA')
legend("bottomright", legend=c("1980-2010", "2030-60", "2070-99"), col=c(wesgreen,wesblue,wesred),lty=1)


dev.off()
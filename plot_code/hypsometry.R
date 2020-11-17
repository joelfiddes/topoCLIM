require(raster)
glacier =shapefile("/home/joel/data/glims/alps/rgi/11_rgi60_CentralEurope.shp")
ele = raster("/home/joel/sim/qmap/ch_tmapp2/eleCH.tif")
sweHIST = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/HIST_2__map.tif")
sweHIST[sweHIST>3]<-3

swercp26near = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP26_2_2030-01-01_2050-12-31__map.tif")
swercp26far = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP26_2_2080-01-01_2099-12-31__map.tif")


swercp85near = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP85_2_2030-01-01_2050-12-31__map.tif")
swercp85far = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP85_2_2080-01-01_2099-12-31__map.tif")






#=====================================================================================
# NEW CODE
#=====================================================================================

ele = raster("/home/joel/sim/qmap/ch_tmapp2/eleCH.tif")



#=====================================================================================
# Hist
#=====================================================================================
step=100
llvec=seq(400,4800,step)
el=llvec[1:length(llvec)-1]

swe_raster = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/HIST_2__map.tif")
swe_raster[swe_raster>3]<-3

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
swe_raster = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP85_2_2030-01-01_2050-12-31__map.tif")
swe_raster[swe_raster>3]<-3

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
swe_raster = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP85_2_2080-01-01_2099-12-31__map.tif")
swe_raster[swe_raster>3]<-3

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
swe_raster = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP26_2_2030-01-01_2050-12-31__map.tif")
swe_raster[swe_raster>3]<-3

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
swe_raster = raster("/home/joel/sim/qmap/ch_tmapp2/ONEMODEL_RESULTS/RCP26_2_2080-01-01_2099-12-31__map.tif")
swe_raster[swe_raster>3]<-3

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


par(mfrow=c(1,2))
lwd=2


#rcp26
mymean = ml1
mysd = sl1
plot(mymean,el,lwd=lwd, xlim=c(0,3), typ="l", col="green")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=rgb(0, 1, 0,0.5),border='NA')

mymean = ml26nr
mysd = sl26nr
lines(mymean,el,lwd=lwd,col="blue")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=rgb(0, 0, 1,0.5),border='NA')


mymean = ml26fr
mysd = sl26fr
lines(mymean,el,lwd=lwd,col="red")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=rgb(1, 0, 0,0.5),border='NA')

# rcp 85
mymean = ml1
mysd = sl1
plot(mymean,el,lwd=lwd, xlim=c(0,3), typ="l", col="green")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=rgb(0, 1, 0,0.5),border='NA')

mymean = ml85nr
mysd = sl85nr
lines(mymean,el,lwd=lwd, col="blue")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=rgb(0, 0, 1,0.5),border='NA')


mymean = ml85fr
mysd = sl85fr
lines(mymean,el,lwd=lwd, col="red")
x = c(mymean- mysd,rev(mymean+ mysd))
y=c((el),rev(el))
polygon (x,y, col=rgb(1, 0, 0,0.5),border='NA')




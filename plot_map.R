#R code
args = commandArgs(trailingOnly=TRUE)
root = args[1]
landformfile = args[2]
maxhs = as.numeric(args[3])

#crispSpatialNow<-function(resultsVec, landform){
pdf(paste0(root,"spatial_plots.pdf"))
par(mfrow=c(2,3))
require(viridis)
require(raster)
#root="/home/joel/sim/qmap/topoclim_test_hpc/g3"
zlim=c(0,1.5)	
mycol=viridis(100)
landform=raster(landformfile)


		file = "meanhist.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		hist=round(rst,2)
		plot(hist,zlim=zlim, col=mycol, main= "1979-2006 Hist ")
	
		file = "mean2030.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot(rst,zlim=zlim, col=mycol, main= "2030-50 RCP26 ")

		file = "mean2080.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)

		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot(rst, zlim=zlim, col=mycol,main= "2080-99 RCP26 ")


		file = "mean2030_rcp85.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot(rst, zlim=zlim, col=mycol,main= "2030-40 RCP85 ")

		file = "mean2080_rcp85.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot(rst, zlim=zlim, col=mycol,main= "2080-99 RCP85 ")

		dev.off()
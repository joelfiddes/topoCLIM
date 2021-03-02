#R code
args = commandArgs(trailingOnly=TRUE)
root = args[1]
landformfile = args[2]
maxhs = as.numeric(args[3])

#crispSpatialNow<-function(resultsVec, landform){
pdf(paste0(root,"/spatial_plots.pdf"))
par(mfrow=c(2,3))
require(viridis)
require(raster)
#root="/home/joel/sim/qmap/topoclim_test_hpc/g3"
zlim=c(0,maxhs)	
mycol=viridis(100)
landform=raster(landformfile)


		file = "meanera5.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		hist=round(rst,2)
		plot(hist,zlim=zlim, col=mycol, main= "1980-2000 Hist ")
		writeRaster(hist, paste0(root,"/era5_spatial.tif"),overwrite=TRUE)

		dev.off()
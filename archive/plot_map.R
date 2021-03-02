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
zlim=c(-1,1)
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
		plot(hist,zlim=zlim, col=mycol, main= "1980-2000 Hist ")
		writeRaster(hist, paste0(root,"/hist.tif"),overwrite=TRUE)
	
		file = "mean2030.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot((rst-hist)/hist,zlim=zlim, col=mycol, main= "2030-50 RCP26 ")
		writeRaster(rst, paste0(root,"/rcp26_2030.tif"),overwrite=TRUE)

		file = "mean2080.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)

		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot((rst-hist)/hist, zlim=zlim, col=mycol,main= "2080-99 RCP26 ")
		writeRaster(rst, paste0(root,"/rcp26_2080.tif"),overwrite=TRUE)

		file = "meanhist.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)

		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		hist=round(rst,2)
		plot(hist,zlim=zlim, col=mycol, main= "1980-2000 Hist ")

		file = "mean2030_rcp85.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot((rst-hist)/hist, zlim=zlim, col=mycol,main= "2030-40 RCP85 ")
		writeRaster(rst, paste0(root,"/rcp85_2030.tif"),overwrite=TRUE)

		file = "mean2080_rcp85.csv"
		require(raster)
		infile =read.csv(paste0(root, "/", file), header=F)
		
		resultsVec=infile$V2  
		l <- length(resultsVec)
		s <- 1:l
		df <- data.frame(s,resultsVec)
		rst <- subs(landform, df,by=1, which=2)
		rst=round(rst,2)
		plot((rst-hist)/hist, zlim=zlim, col=mycol,main= "2080-99 RCP85 ")
		writeRaster(rst, paste0(root,"/rcp85_2080.tif"),overwrite=TRUE)

		dev.off()



		# require(raster)
		# rasters1 <- list.files(".", pattern="*.tif",full.names=TRUE, recursive=FALSE)
		# rast.list <- list()
		#   for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

		# # And then use do.call on the list of raster objects
		# rast.list$fun <- mean
		# rast.mosaic <- do.call(mosaic,rast.list)
		# rst<-rast.mosaic
		# writeRaster(rst, "hansenGlobal.tif", overwrite=T)
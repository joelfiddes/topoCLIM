require(raster)

args = commandArgs(trailingOnly=TRUE)
tsub_root = args[1]
meanVar= args[2]
nclust = as.numeric(args[3])
outname = args[4]

gridseq =c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,3,4,5,6,7,8,9)

meanswe =read.csv(meanVar ,head=F )


# this is given by filenames = sorted(glob.glob(wd + "/*/*/listpoints.txt"))    
# in tscale3D.py but only valid for this domain - NEED TO FIX!

	for (i in 1:length(gridseq)){
		print(i)
		grid = gridseq[i]
		sample_indexes = (i-1)*nclust+(1:nclust)

		swe =meanswe[sample_indexes,2]
		s <- 1:length(swe)
		df <- data.frame(s,swe)
		landform=raster(paste0(tsub_root,"/sim/g",grid,"/landform.tif"))

		rst <- subs(landform, df,by=1, which=2)
		writeRaster(rst, paste0(tsub_root, "/g",grid,"_swe.tif"),overwrite=TRUE)
		}





require(raster)
require(viridis)
	
zlim=c(0,3)	
mycol=viridis(100)

pdf(paste0(tsub_root,"/",outname,"_map.pdf"))

rasters1 <- list.files(tsub_root, pattern="*swe.tif",full.names=TRUE, recursive=TRUE)
rast.list <- list()
for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

	# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
hist<-rast.mosaic
	#plot(hist,zlim=zlim, col=mycol, main= "1980-2000 Hist ")
writeRaster(hist, paste0(tsub_root,"/",outname,"__map.tif"), overwrite=T)
plot(hist, zlim=zlim, col=mycol)
dev.off()
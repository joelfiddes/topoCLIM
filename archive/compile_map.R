args = commandArgs(trailingOnly=TRUE)
root = args[1]
outname=args[2]
maxhs=as.numeric(args[3])

require(raster)
require(viridis)
zlim=c(0,maxhs)	
zlim=c(-1,1)	
mycol=viridis(10)

pdf(paste0(root,"/",outname,"_spatial_plots_compiled.pdf"))
par(mfrow=c(2,2))
rasters1 <- list.files(root, pattern="*hist.tif",full.names=TRUE, recursive=TRUE)

rast.list <- list()
  for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
hist<-rast.mosaic
#plot(hist,zlim=zlim, col=mycol, main= "1980-2000 Hist ")
writeRaster(hist, paste0(root,"/",outname,"_hist_comp.tif"), overwrite=T)


rasters1 <- list.files(root, pattern="*rcp26_2030.tif",full.names=TRUE, recursive=TRUE)

rast.list <- list()
  for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
rst<-rast.mosaic
plot((rst-hist)/hist,zlim=zlim, col=mycol, main= "2030-50 RCP26 ")
writeRaster(rst, paste0(root,"/",outname,"_rcp26_2030_comp.tif"), overwrite=T)

rasters1 <- list.files(root, pattern="*rcp26_2080.tif",full.names=TRUE, recursive=TRUE)

rast.list <- list()
  for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
rst<-rast.mosaic

		plot(rst, zlim=zlim, col=mycol,main= "2080-99 RCP26 ")
writeRaster((rst-hist)/hist, paste0(root,"/",outname,"_rcp26_2080_comp.tif"), overwrite=T)

rasters1 <- list.files(root, pattern="*rcp85_2030.tif",full.names=TRUE, recursive=TRUE)

rast.list <- list()
  for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
rst<-rast.mosaic
plot((rst-hist)/hist, zlim=zlim, col=mycol,main= "2030-40 RCP85 ")
writeRaster(rst, paste0(root,"/",outname,"_rcp85_2030_comp.tif"), overwrite=T)


rasters1 <- list.files(root, pattern="*rcp85_2080.tif",full.names=TRUE, recursive=TRUE)

rast.list <- list()
  for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
rst<-rast.mosaic
plot((rst-hist)/hist, zlim=zlim, col=mycol,main= "2080-99 RCP85 ")
writeRaster(rst, paste0(root,"/",outname,"_rcp85_2080_comp.tif"), overwrite=T)
dev.off()
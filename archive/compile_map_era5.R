args = commandArgs(trailingOnly=TRUE)
root = args[1]
outname=args[2]
maxhs=as.numeric(args[3])

require(raster)
require(viridis)
zlim=c(0,maxhs)	
mycol=viridis(10)

pdf(paste0(root,"/",outname,"_spatial_plots_compiled.pdf"))

rasters1 <- list.files(root, pattern="era5_spatial.tif",full.names=TRUE, recursive=TRUE)

rast.list <- list()
  for(i in 1:length(rasters1)) { rast.list[i] <- raster(rasters1[i]) }

# And then use do.call on the list of raster objects
rast.list$fun <- mean
rast.mosaic <- do.call(mosaic,rast.list)
rst<-rast.mosaic
plot(rst,zlim=zlim, col=mycol, main= "1980-2000 Hist ")
writeRaster(rst, paste0(root,"/",outname,"_hist_comp.tif"), overwrite=T)
dev.off()
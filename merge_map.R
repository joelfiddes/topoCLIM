require(raster)
require(viridis)

dem=raster("/home/joel/sim/qmap/ch_tmapp2/predictors/ele.tif")
glaciers=shapefile("/home/joel/data/glims/alps/glims_download_17927/glims_polygons.shp")
gl = crop(glaciers,ch)
hist_swe = raster("/home/joel/sim/qmap/ch_tmapp2/HIST_2__map.tif")
rcp26near_swe = raster("/home/joel/sim/qmap/ch_tmapp2/RCP26_2_2030-01-01_2050-12-31__map.tif")
rcp26far_swe = raster("/home/joel/sim/qmap/ch_tmapp2/RCP26_2_2080-01-01_2099-12-31__map.tif")
rcp85near_swe = raster("/home/joel/sim/qmap/ch_tmapp2/RCP85_2_2030-01-01_2050-12-31__map.tif")
rcp85far_swe = raster("/home/joel/sim/qmap/ch_tmapp2/RCP85_2_2080-01-01_2099-12-31__map.tif")

hist_swe[hist_swe==0]<-NA
rcp26near_swe[rcp26near_swe==0]<-NA
rcp26far_swe[rcp26far_swe==0]<-NA
rcp85near_swe[rcp85near_swe==0]<-NA
rcp85far_swe[rcp85far_swe==0]<-NA

maxhs=0.7
hist_swe[hist_swe>maxhs]<-NA
rcp26near_swe[rcp26near_swe>maxhs]<-NA
rcp26far_swe[rcp26far_swe>maxhs]<-NA
rcp85near_swe[rcp85near_swe>maxhs]<-NA
rcp85far_swe[rcp85far_swe>maxhs]<-NA

demch = mask(crop(dem,ch),ch)
hist_swe_ch = mask(crop(hist_swe,ch),ch)
rcp26near_swe_ch = mask(crop(rcp26near_swe,ch),ch)
rcp26far_swe_ch = mask(crop(rcp26far_swe,ch),ch)
rcp85near_swe_ch = mask(crop(rcp85near_swe,ch),ch)
rcp85far_swe_ch = mask(crop(rcp85far_swe,ch),ch)


# zlim=c(-5,0)	
# mycol=inferno(100)
# par(mfrow=c(2,3))
# plot(demch, legend=F, col=mygrey(50))
# plot(hist_swe_ch, col=mycol, zlim=zlim, add=T)
# plot(demch, legend=F, col=mygrey(50))
# plot(rcp26near_swe_ch, col=mycol, zlim=zlim,add=T)
# plot(demch, legend=F, col=mygrey(50))
# plot(rcp26far_swe_ch, col=mycol, zlim=zlim,add=T)
# plot(demch, legend=F, col=mygrey(50))
# plot(hist_swe_ch, col=mycol, zlim=zlim,add=T)
# plot(demch, legend=F, col=mygrey(50))
# plot(rcp85near_swe_ch, col=mycol, zlim=zlim,add=T)
# plot(demch, legend=F, col=mygrey(50))
# plot(rcp85far_swe_ch, col=mycol, zlim=zlim,add=T)
library(colorspace)
myTheme <- rasterTheme(region=viridis(10))
demstk=stack(demch, demch, demch,demch, demch, demch)
stk = stack(hist_swe_ch, hist_swe_ch, rcp26near_swe_ch,  rcp85near_swe_ch,rcp26far_swe_ch, rcp85far_swe_ch)
names(demstk)<- c("Hist_1980_00", "Hist_1980_2000", "RCP26_2030_50","RCP85_2030_50", "RCP26_2080_99", "RCP85_2080_99")
names(stk)<- c("Hist_1980_00", "Hist_1980_2000", "RCP26_2030_50","RCP85_2030_50", "RCP26_2080_99", "RCP85_2080_99")
mycol=viridis(100)
p0 = levelplot(demstk, maxpixels=1e6,par.settings = GrTheme)
p1 = levelplot(stk, par.settings = myTheme)



png("/home/joel/manuscripts/qmap/plots/swe.png", height=600, width=800)
p1 + as.layer(p0, under = TRUE) #+ layer(sp.lines(gl, lwd=0.2, col='red'))
dev.off()





zlim=c(0,1)	
mycol=viridis(100)
par(mfrow=c(2,3))
plot(hist, col=mycol, zlim=zlim)
plot(glaciers, col='red', add=T)
plot(rcp26near, col=mycol, zlim=zlim)
plot(rcp26far, col=mycol, zlim=zlim)
plot(hist, col=mycol, zlim=zlim)
plot(rcp85near, col=mycol, zlim=zlim)
plot(rcp85far, col=mycol, zlim=zlim)


zlim=c(-1,0)	
mycol=viridis(100)
par(mfrow=c(2,2))
plot((rcp26near-hist), col=mycol, zlim=zlim)
plot(rcp26far-hist, col=mycol, zlim=zlim)
plot(rcp85near-hist, col=mycol, zlim=zlim)
plot(rcp85far-hist, col=mycol, zlim=zlim)


#=========================================================================================
# PERMAFROST
#=========================================================================================

ch =getData('GADM', country='CHE', level=0)  
mygrey <- colorRampPalette(c("grey0", "grey100"))

require(raster)
require(viridis)

shp = 
hist_gst = raster("/home/joel/sim/qmap/ch_tmapp2/HIST_5__map.tif")
rcp26near_gst = raster("/home/joel/sim/qmap/ch_tmapp2/RCP26_5_2030-01-01_2050-12-31__map.tif")
rcp26far_gst = raster("/home/joel/sim/qmap/ch_tmapp2/RCP26_5_2080-01-01_2099-12-31__map.tif")
rcp85near_gst = raster("/home/joel/sim/qmap/ch_tmapp2/RCP85_5_2030-01-01_2050-12-31__map.tif")
rcp85far_gst = raster("/home/joel/sim/qmap/ch_tmapp2/RCP85_5_2080-01-01_2099-12-31__map.tif")

hist_gst[hist_gst>0]<-NA
rcp26near_gst[rcp26near_gst>0]<-NA
rcp26far_gst[rcp26far_gst>0]<-NA
rcp85near_gst[rcp85near_gst>0]<-NA
rcp85far_gst[rcp85far_gst>0]<-NA


demch = mask(crop(dem,ch),ch)
hist_gst_ch = mask(crop(hist_gst,ch),ch)
rcp26near_gst_ch = mask(crop(rcp26near_gst,ch),ch)
rcp26far_gst_ch = mask(crop(rcp26far_gst,ch),ch)
rcp85near_gst_ch = mask(crop(rcp85near_gst,ch),ch)
rcp85far_gst_ch = mask(crop(rcp85far_gst,ch),ch)
#glac_ch  = mask(crop(glaciers,ch),ch)


zlim=c(-5,0)	
mycol=inferno(100)
par(mfrow=c(2,3))
plot(demch, legend=F, col=mygrey(50))
plot(hist_gst_ch, col=mycol, zlim=zlim, add=T)
plot(demch, legend=F, col=mygrey(50))
plot(rcp26near_gst_ch, col=mycol, zlim=zlim,add=T)
plot(demch, legend=F, col=mygrey(50))
plot(rcp26far_gst_ch, col=mycol, zlim=zlim,add=T)
plot(demch, legend=F, col=mygrey(50))
plot(hist_gst_ch, col=mycol, zlim=zlim,add=T)
plot(demch, legend=F, col=mygrey(50))
plot(rcp85near_gst_ch, col=mycol, zlim=zlim,add=T)
plot(demch, legend=F, col=mygrey(50))
plot(rcp85far_gst_ch, col=mycol, zlim=zlim,add=T)

demstk=stack(demch, demch, demch,demch, demch, demch)
stk = stack(hist_gst_ch,hist_gst_ch,  rcp26near_gst_ch, rcp85near_gst_ch,rcp26far_gst_ch, rcp85far_gst_ch)

names(demstk)<- c("Hist_1980_00", "Hist_1980_2000", "RCP26_2030_50","RCP85_2030_50", "RCP26_2080_99", "RCP85_2080_99")
names(stk)<- c("Hist_1980_00", "Hist_1980_2000", "RCP26_2030_50","RCP85_2030_50", "RCP26_2080_99", "RCP85_2080_99")
p0 = levelplot(demstk, maxpixels=1e6,par.settings = GrTheme)
p1 = levelplot(stk, zlim=zlim)

png("/home/joel/manuscripts/qmap/plots/GST.png", height=600, width=800)
p1 + as.layer(p0, under = TRUE) #+ layer(sp.lines(glaciers, lwd=0.2, col='lightblue'))
dev.off()
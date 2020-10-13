require(raster)
r=raster("/home/joel/sim/qmap/GR_data/forcing/SURF_timmean.nc", var='tp')
snow=raster('/home/joel/sim/qmap/topoclim_test_hpc/swe_hist_comp.tif')

l2 = crop(r,snow)


disagg = res(r)[1]/res(snow)[1]

l1 = disaggregate(r, disagg)
l = disaggregate(r, disagg, method='bilinear')

snowscale = l/l1

snowscale_resamp = resample(snowscale, snow)


par(mfrow=c(1,2))                                                                                                                                                                                                                                 
plot(snowscale_resamp*snow)
plot(snow)  
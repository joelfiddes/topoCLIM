#constructs api requests and renames wget files to keep tidy
# variable


library(RCurl)

# wind
# # fixed:
# domain='EUR-44'
# freq='6hr' # '3hr' or 'day'
# outpath='/home/joel/sim/qmap/'
# # iterate over:
#vars=c( 'uas', 'vas')
# expers=c('historical', 'rcp26', 'rcp85') #
#
# for (var in vars){
#   for (exper in expers){
#   APIpattern=
#   paste0("https://esgf-data.dkrz.de/esg-search/wget?project=CORDEX&variable=",var,"&time_frequency=",freq,"&domain=",domain,"&experiment=",exper,"&download_structure=project,product,domain,institute,driving_model,experiment,ensemble,rcm_name,rcm_version,time_frequency,variable")
#   print(var)
#   print(exper)
#   write(getURL(APIpattern),paste0(outpath,domain,"_",freq,"_", var,"_",exper,".sh"))
#   }
#
# }



# fixed:
domain='EUR-44'
freq='day' # '3hr' or 'day'
outpath='/home/joel/sim/qmap/test/'
indir=paste0(outpath,"CORDEX/output/", domain,"/")
res = 0.5 # output grid resolution in degrees
setwd(outpath)
# iterate over:
vars=c('tas', 'pr', 'ps', 'hurs', 'rsds','rlds', 'uas', 'vas')
expers=c('historical', 'rcp26', 'rcp85') #

# postprocessing
# variable



# define domain here:
lonE <- 5
lonW<-11
latS<-45
latN<-48

# define and construct coords config for remapbil
coordsPath <- paste0(outpath, "/coords.txt") # describes reprojection

# params to write
gridtype <- 'lonlat'
xsize <- (lonW-lonE)/res
ysize <- (latN-latS)/res
xfirst <- lonE
xinc <- res
yfirst <- latS
yinc <- res

# write remapbil config file
write(paste0("gridtype = ",gridtype),paste0(outpath, "/coords.txt"))
write(paste0("xsize = ",xsize),paste0(outpath, "/coords.txt"), append=T)
write(paste0("ysize = ",ysize),paste0(outpath, "/coords.txt"), append=T)
write(paste0("xfirst = ",xfirst),paste0(outpath, "/coords.txt"), append=T)
write(paste0("xinc = ",xinc),paste0(outpath, "/coords.txt"), append=T)
write(paste0("yfirst = ",yfirst),paste0(outpath, "/coords.txt"), append=T)
write(paste0("yinc = ",yinc),paste0(outpath, "/coords.txt"), append=T)

# run code!

for (exper in expers){
  for (var in vars){
  APIpattern=
  paste0("https://esgf-data.dkrz.de/esg-search/wget?project=CORDEX&institute=CNRM","&variable=",var,"&time_frequency=",freq,"&domain=",domain,"&experiment=",exper,"&download_structure=project,product,domain,institute,driving_model,experiment,ensemble,rcm_name,rcm_version,time_frequency,variable")
  print(APIpattern)
  write(getURL(APIpattern),paste0(outpath,domain,"_",freq,"_", var,"_",exper,".sh"))

  # chmod
  system(paste0("chmod u+x ", outpath,domain,"_",freq,"_", var,"_",exper,".sh"))

  # run download
  system(paste0("bash ", outpath,domain,"_",freq,"_", var,"_",exper,".sh"))

  # find all files to postprocess
  filesAll <- list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)
  freq_index <- grep(freq, filesAll)
  files <- filesAll[freq_index]
  donefiles = list.files(path=indir, pattern="ll.nc$", recursive=T, full.name=T)
  TSfiles = list.files(path=indir, pattern="TS.nc$", recursive=T, full.name=T)
  donefilesIndex = which(!files %in% donefiles)

  filestodo = files[donefilesIndex]
  tsindex = which(filestodo %in% TSfiles)

# check if there any TSfiles, in vector length is ) (none) then this would give no results with a -0 index
  if (length(tsindex > 0)){
    filestoreallydo <- filestodo[-tsindex]
    }
  if ((length(tsindex) == 0)){
    filestoreallydo <- filestodo
    }
  # exclude all _ll.nc files here in case of reprocessing otherwise get .._ll_ll.nc files!
  # exclude all _TS.nc files here otherwise removes .._TS.nc files and makes TS_ll.nc files!
    for ( file in filestoreallydo ) {

    	basename= strsplit(file,".nc")
    	#system(paste0("cdo remapbil,",coordsPath," -sellonlatbox,",lonE,",",lonW,",",latS,",",latN," ", file," " ,basename,"_ll.nc"))
    	system(paste0("cdo remapbil,",coordsPath," ", file," " ,basename,"_ll.nc")) # just need to remapbil to smaller domain and implicitly cuts
    	system(paste0("rm ", file))
    }


  }

}

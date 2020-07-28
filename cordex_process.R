# source
#https://stackoverflow.com/questions/31549880/how-to-convert-a-rotated-netcdf-back-to-normal-lat-long

# recursively searchs indir for nc files
# cuts to bbox
# remaps all rotated grids to ll according coordsPath config
# removes all original files
# concatentes individual periods to single file _TS.nc
# requires CDO

# variable
indir="/home/joel/sim/qmap/CORDEX/output/EUR-44/"
res = 0.5 # output grid resolution in degrees
freq='3hr' #'day' or '3hr' to subset day my temporal frequency

# define domain here:
lonE <- 5
lonW<-11
latS<-45
latN<-48

# define and construct coords config for remapbil
coordsPath <- paste0(indir, "/coords.txt") # describes reprojection

# params to write
gridtype <- 'lonlat'
xsize <- (lonW-lonE)/res
ysize <- (latN-latS)/res
xfirst <- lonE
xinc <- res
yfirst <- latS
yinc <- res


# write remapbil config file
write(paste0("gridtype = ",gridtype),paste0(indir, "/coords.txt"))
write(paste0("xsize = ",xsize),paste0(indir, "/coords.txt"), append=T)
write(paste0("ysize = ",ysize),paste0(indir, "/coords.txt"), append=T)
write(paste0("xfirst = ",xfirst),paste0(indir, "/coords.txt"), append=T)
write(paste0("xinc = ",xinc),paste0(indir, "/coords.txt"), append=T)
write(paste0("yfirst = ",yfirst),paste0(indir, "/coords.txt"), append=T)
write(paste0("yinc = ",yinc),paste0(indir, "/coords.txt"), append=T)

#code
filesAll <- list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)
freq_index <- grep(freq, filesAll)
files <- filesAll[freq_index]
donefiles = list.files(path=indir, pattern="ll.nc$", recursive=T, full.name=T)
donefilesIndex = which(!files %in% donefiles)
filestodo = files[donefilesIndex]
# exclude all _ll.nc files here in case of reprocessing otherwise get .._ll_ll.nc files!
for ( file in filestodo ) {

	basename= strsplit(file,".nc")
	#system(paste0("cdo remapbil,",coordsPath," -sellonlatbox,",lonE,",",lonW,",",latS,",",latN," ", file," " ,basename,"_ll.nc"))
	system(paste0("cdo remapbil,",coordsPath," ", file," " ,basename,"_ll.nc")) # just need to remapbil to smaller domain and implicitly cuts
	#system(paste0("rm ", file))
}

# delete all original downloads 100GB ->
filesProcessed = list.files(path=indir, pattern="_ll.nc$", recursive=T, full.name=T)
filesAll = list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)
keepIndex = which(filesAll %in% filesProcessed)
files2keep = filesAll[keepIndex]

rmIndex = which(!filesAll %in% filesProcessed)
files2rm = filesAll[rmIndex]

for (i in files2rm){
	system(paste0("rm ", i))
}

if(freq=='day'){strsplitPar <- 23}
if(freq=='3hr'){strsplitPar <- 31}
files = list.files(path=indir, pattern="_ll.nc$", recursive=T, full.name=T)

base = substring(files, 1,nchar(files)-strsplitPar)
base2 = unique(base)
for (mybase in base2){
#mybase2 = unlist(strsplit(mybase,'/'))[[length(unlist(strsplit(mybase,'/')))]]
#base3 = list.files(path=indir, pattern=mybase2, recursive=T, full.name=T)
system(paste0("cdo -b F64 -f nc2 mergetime ",mybase,"* ", mybase,"TS.nc"))
print(paste0("concatenated file: ", mybase,"TS.nc"))
}

final_TS = list.files(path=indir, pattern="_TS.nc$", recursive=T, full.name=T)

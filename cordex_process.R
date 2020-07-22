# source
#https://stackoverflow.com/questions/31549880/how-to-convert-a-rotated-netcdf-back-to-normal-lat-long

# recursively searchs indir for nc files
# cuts to bbox
# remaps all rotated grids to ll according coordsPath config
# removes all original files
# concatentes individual periods to single file _TS.nc
# requires CDO

# variable
res=0.5
indir="/home/joel/sim/qmap/CORDEX/output/EUR-44/"
coordsPath = paste0(indir, "/coords.txt") # describes reprojection
lonE = 5
lonW=11
latS=45
latN=48

# construct coords config
gridtype = 'lonlat'
xsize = (lonW-lonE)/res
ysize = (latN-latS)/res
xfirst = lonE
xinc = res
yfirst = latS
yinc = res


write(paste0("gridtype = ",gridtype),paste0(indir, "/coords.txt"))
write(paste0("xsize = ",xsize),paste0(indir, "/coords.txt"), append=T)
write(paste0("ysize = ",ysize),paste0(indir, "/coords.txt"), append=T)
write(paste0("xfirst = ",xfirst),paste0(indir, "/coords.txt"), append=T)
write(paste0("xinc = ",xinc),paste0(indir, "/coords.txt"), append=T)
write(paste0("yfirst = ",yfirst),paste0(indir, "/coords.txt"), append=T)
write(paste0("yinc = ",yinc),paste0(indir, "/coords.txt"), append=T)

file="/home/joel/sim/qmap/CORDEX/output/EUR-44/KNMI/ICHEC-EC-EARTH/rcp26/r12i1p1/RACMO22E/v1/day/pr/pr_EUR-44_ICHEC-EC-EARTH_rcp26_r12i1p1_KNMI-RACMO22E_v1_day_20060101-20101231.nc"
	basename= strsplit(file,".nc")
system(paste0("cdo sellonlatbox,",lonE,",",lonW,",",latS,",",latN," ", file," " ,basename,"_SUB.nc"))
system(paste0("cdo remapbil,",coordsPath," ", file," " ,basename,"_ll3.nc"))

system(paste0("cdo remapbil,",coordsPath," -sellonlatbox,",lonE,",",lonW,",",latS,",",latN," ", file," " ,basename,"_ll.nc"))
system(paste0("cdo sellonlatbox,",lonE,",",lonW,",",latS,",",latN," -remapbil,",coordsPath," ", file," " ,basename,"_ll2.nc"))
#code
files = list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)

# exclude all _ll.nc files here in case of reprocessing otherwise get .._ll_ll.nc files!
for ( file in files ) {

	basename= strsplit(file,".nc")
	#system(paste0("cdo remapbil,",coordsPath," -sellonlatbox,",lonE,",",lonW,",",latS,",",latN," ", file," " ,basename,"_ll.nc"))
	#system(paste0("cdo sellonlatbox,",lonE,",",lonW,",",latS,",",latN," -remapbil,",coordsPath," ", file," " ,basename,"_ll2.nc")) #we need to reproject and THEN cut to get correct results (takes 0.04s longer - not signif)
	system(paste0("cdo remapbil,",coordsPath," ", file," " ,basename,"_ll.nc")) # just need to remapbil to smaller domain and implicitly cuts
	system(paste0("rm ", file))
}


files = list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)

base = substring(files, 1,nchar(files)-23)
base2 = unique(base)
for (mybase in base2){
#mybase2 = unlist(strsplit(mybase,'/'))[[length(unlist(strsplit(mybase,'/')))]]
#base3 = list.files(path=indir, pattern=mybase2, recursive=T, full.name=T)
system(paste0("cdo -b F64 -f nc2 mergetime ",mybase,"* ", mybase,"TS.nc"))
print(paste0("concatenated file: ", mybase,"TS.nc"))
}

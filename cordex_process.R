# source https://stackoverflow.com/questions/31549880/how-to-convert-a-rotated-netcdf-back-to-normal-lat-long

# recursively searchs indir for nc files
# cuts to bbox
# remaps all rotated grids to ll according coordsPath config
# removes all original files
# concatentes individual periods to single file _TS.nc
# requires CDO

# variable
indir="/home/joel/sim/qmap/CORDEX/"
coordsPath = "/home/joel/src/topoCLIM/coords.txt" # describes reprojection
lonE = 5.548096
lonW=10.601807
latS=45.460131
latN=47.945786

#code
files = list.files(path=indir, pattern=".nc$", recursive=T, full.name=T)

dir.create(outdir)

for ( file in files ) {

	basename= strsplit(file,".nc")
	system(paste0("cdo remapbil,",coordsPath," -sellonlatbox,",lonE,",",lonW,",",latS,",",latN," ", file," " ,basename,"_ll.nc"))
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



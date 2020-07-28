import glob
import os
outdir="/home/joel/sim/qmap/test/pyout/"
domain='EUR-44'

# define domain here:
lonE = int(5)
lonW = int(11)
latS = int(45)
latN = int(48)
res=0.5 # output grid resolution in degrees

# define and construct coords config for remapbil which  describes reprojection
coordsPath = outdir+ "coords.txt"

# write remapbil config file
f = open(coordsPath, "w")
f.write("gridtype = " + "lonlat\n")
f.write("xsize = " + str(int((lonW-lonE)/res))+"\n")
f.write("ysize = " + str(int((latN-latS)/res))+"\n")
f.write("xfirst = " + str(lonE)+"\n")
f.write("xinc = " + str(res)+"\n")
f.write("yfirst = " + str(latS)+"\n")
f.write("yinc = "+ str(res)+"\n")
f.close()

ts_dir=outdir+"TS/"
if not os.path.exists(ts_dir):
	os.makedirs(ts_dir)

mydownloads = glob.glob(outdir+"*.nc")
mydownloads.sort()
root  = [i.split('day', 1)[0] for i in mydownloads]
root_unique = list(set(root))

for ru in root_unique:
    basename = os.path.basename(ru)
    os.system("cdo -b F64 -f nc2 mergetime "+ru+"* "+ts_dir+basename+"_TS.nc")
    print(("concatenated file: "+ ts_dir+ru,"_TS.nc"))


# merg on vars
TSfiles = glob.glob(ts_dir+"*_TS.nc")
root  = [i.split(domain+"_", 1)[1] for i in TSfiles]
root_unique = list(set(root))

for ru in root_unique:
    os.system("cdo merge "+ ts_dir+"*"+ru +" "+ts_dir+ru+"_TS_ALL.nc")
    os.system("cdo remapbil," + coordsPath+" "+ ts_dir+ru+"_TS_ALL.nc"+" " +ts_dir+ru+"_TS_ALL_ll.nc")


#cleanup
os.system("rm "+ ts_dir+"*_TS_ALL.nc")
os.system("rm "+ ts_dir+"*_TS.nc")

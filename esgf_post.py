import glob
import os
import xarray as xr
import pandas as pd

outdir="/home/joel/sim/qmap/test/pyout/"
domain='EUR-44'
time_frequency='day' #3hr' #
ts_dir=outdir+"aresult/"

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

# time merge
if not os.path.exists(ts_dir):
	os.makedirs(ts_dir)

mydownloads = glob.glob(outdir+"*.nc")
mydownloads.sort()
root  = [i.split(time_frequency, 1)[0] for i in mydownloads]
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

# analyse results
results = glob.glob(ts_dir+"*.nc")



d=[]

for my_file in results:
	ds = xr.open_dataset(my_file)
	basename = os.path.basename(my_file)
	start = ds['time'][1].data  
	end = ds['time'][-1].data  
	nvars = len(ds.data_vars) 

	d.append(
	{
	'File': basename,
	'Start': start,
	'End': end,
	'Variables': nvars

	}


		)
df =pd.DataFrame(d)[['File', 'Start', 'End', 'Variables']]

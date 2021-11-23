"""
    DESCRIPTION:
        This helper function postprocesses downloaded CORDEX data. (1) 
        reprojects to long lat (2) crop to domain defined by 
        "mydomain" (3) merges data by time (4) merges data by variable

    ARGS:
        cordex_domain (str): cordex experiment domain e.g. "EUR-44"
        outdir (str): Path tpo write results to e.g. /path/to/results
        mydomain (str): Path to shapefile which defines your domain
        output_res  (str): output resolution of grid (degrees)

    RETURNS:
        NULL (files written to outdir)
    


   """

import glob
import os
import xarray as xr
import pandas as pd
import geopandas as gpd

cordex_domain = sys.argv[1] # CORDEX domain eg "EUR-44"
outdir = sys.argv[2] # where to write results eg /path/to/results
mydomain = sys.argv[3] # path to shapefile which defines your domain
output_res = sys.argv[4] # output resolution of grid (degrees)

# fixed
time_frequency = 'day'  # 3hr' #
ts_dir = outdir + "aresult/"

# define domain here based on era5 domain
era5ref = gpd.read_file(mydomain).total_bounds
lonE = era5ref[2]  # int(5)
lonW = era5ref[0]  # int(11)
latS = era5ref[1]  # int(45)
latN = era5ref[3]  # int(48)
res = output_res # output grid resolution in degrees

# define and construct coords config for remapbil which  describes reprojection
coordsPath = outdir + "coords.txt"

# write remapbil config file
# remapping is done to grid centres thats why we add half a gridbox (res/2)


# If the poly is good this works
f = open(coordsPath, "w")
f.write("gridtype = " + "lonlat\n")
f.write("xsize = " + str(int((lonE - lonW) / res)) + "\n")
f.write("ysize = " + str(int((latN - latS) / res)) + "\n")
f.write("xfirst = " + str(lonW + res / 2) + "\n")
f.write("xinc = " + str(res) + "\n")
f.write("yfirst = " + str(latS + res / 2) + "\n")
f.write("yinc = " + str(res) + "\n")
f.close()


# time merge
if not os.path.exists(ts_dir):
    os.makedirs(ts_dir)

mydownloads = sorted(glob.glob(outdir + "*.nc"))
root = [i.split(time_frequency, 1)[0] for i in mydownloads]
root_unique = list(set(root))

for ru in root_unique:
    basename = os.path.basename(ru)
    os.system(
        "cdo -b F64 -f nc2 mergetime " +
        ru +
        "* " +
        ts_dir +
        basename +
        "_TS.nc")
    print(("concatenated file: " + ts_dir + ru, "_TS.nc"))


# merg on vars
TSfiles = glob.glob(ts_dir + "*_TS.nc")
root = [i.split(domain + "_", 1)[1] for i in TSfiles]
root_unique = list(set(root))

for ru in root_unique:
    os.system(
        "cdo merge " +
        ts_dir +
        "*" +
        ru +
        " " +
        ts_dir +
        ru +
        "_TS_ALL.nc")
    os.system(
        "cdo remapbil," +
        coordsPath +
        " " +
        ts_dir +
        ru +
        "_TS_ALL.nc" +
        " " +
        ts_dir +
        ru +
        "_TS_ALL_ll.nc")


# cleanup
os.system("rm " + ts_dir + "*_TS_ALL.nc")
os.system("rm " + ts_dir + "*_TS.nc")


# analyse results
results = glob.glob(ts_dir + "*.nc")

d = []

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
df = pd.DataFrame(d)[['File', 'Start', 'End', 'Variables']]

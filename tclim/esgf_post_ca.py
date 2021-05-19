import glob
import os
import xarray as xr
import pandas as pd
import geopandas as gpd

outdir = "/home/joel/data/cordex/CAS22/"
domain = 'CAS-22'

time_frequency = 'day'  # 3hr' #
ts_dir = outdir + "aresult/"
era5spatialRef = "/home/joel/data/cordex/CAS22_CACdomain.shp"
# define domain here based on era5 domain
era5ref = gpd.read_file(era5spatialRef).total_bounds
lonE = era5ref[2]  # int(5)
lonW = era5ref[0]  # int(11)
latS = era5ref[1]  # int(45)
latN = era5ref[3]  # int(48)
res = 0.22  # output grid resolution in degrees
mydownloads = glob.glob(outdir + "*.nc")

# get long lat
ds = xr.open_dataset(mydownloads[0])
xfirst = ds['lon'][0, :].data[0]  # westmost
yfirst = ds['lat'][:, 0].data[0]  # southernmost
yres = ds['lat'][:, 0].data[1] - yfirst
xres = ds['lon'][0, :].data[1] - xfirst

# define and construct coords config for remapbil which  describes reprojection
coordsPath = outdir + "coords.txt"

# write remapbil config file
# remapping is done to grid centres thats why we add half a gridbox (res/2)


# for cas we do like this, prescribe values due to strange format
# you can get the info by running :
# (base) joel@mountainsense:~$ cdo griddes /home/joel/data/cordex/CAS22/tas_CAS-22_CNRM-CERFACS-CNRM-CM5_rcp45_r1i1p1_RMIB-UGent-ALARO-0_v1_day_20260101-20301231.nc
xsize = len(ds['lon'][0, :].data)
ysize = len(ds['lat'][:, 0].data)

# do we nedd to remap to centres *(ie add half a box to yfirst and xfirst)
f = open(coordsPath, "w")
f.write("gridtype = " + "lonlat\n")
f.write("xsize = " + str(int(xsize)) + "\n")
f.write("ysize = " + str(int(ysize)) + "\n")
f.write("xfirst = " + str(xfirst) + "\n")
f.write("xinc = " + str(xres) + "\n")
f.write("yfirst = " + str(yfirst) + "\n")
f.write("yinc = " + str(yres) + "\n")
f.close()

# time merge
if not os.path.exists(ts_dir):
    os.makedirs(ts_dir)


mydownloads.sort()
root = [i.split(time_frequency, 1)[0] for i in mydownloads]
root_unique = list(set(root))

# cat all file to time


# convert all files to long lat
# os.system(
#         "cdo remapbil," +
#         coordsPath +
#         " " +
#         file +
#         " " +
#         outfile)

# outfile = ts_dir + ru + "_TS_ALL_ll.nc"
# if os.path.isfile(outfile):
#     print(outfile + " already made!")
#     continue

for ru in root_unique:

    basename = os.path.basename(ru)

    outfile = ts_dir + basename + "_TS.nc"
    if os.path.isfile(outfile):
        print(outfile + " already made!")
        os.system("rm " + ru + "* ")
        continue

    if not os.path.isfile(outfile):
        print("not yet!")

        # cdo only works with systema and not subprocess (recommended) for some reason  # noqa: E101
        os.system(
            "cdo -b F64 -f nc2 mergetime " +
            ru +
            "* " +
            ts_dir + basename +
            "_TS.nc")
        print(("concatenated file: " + ts_dir + basename + "_TS.nc"))

    #print("moving: " )
    #print(glob.glob(ru +"*") )
    #os.system("mv " + ru +"* " +" /home/joel/mnt/myserver/nas/data/CORDEX/CAS22")


# merg on vars
TSfiles = glob.glob(ts_dir + "*_TS.nc")

for tfile in TSfiles:
    print(tfile)

    os.system(
        "cdo remapbil," +
        coordsPath +
        " " +
        tfile +
        " " +
        tfile +
        "_ll.nc")

    os.system("rm " + tfile)


# cleanup
#os.system("rm " + ts_dir + "*_TS_ALL.nc")
#os.system("rm " + ts_dir + "*_TS.nc")


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

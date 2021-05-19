import pandas as pd
import glob
import resample_timeseries as resamp
import topoCLIM as tclim
import subprocess
import os
# asuume run:
# esgf_get.py
# esgf_post.py


# ===============================================================================
# topoCLIM.py settings
# ===============================================================================
wd = "/home/joel/sim/qmap/topoclim_test/"
raw_dir = '/home/joel/sim/qmap/test/pyout/'
# standard calender examples (can be any file with standard cal) to interp
# non-standard cals to in hist and clim period
nc_standard_clim = '/home/joel/sim/qmap/test/pyout/aresult/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
nc_standard_hist = '/home/joel/sim/qmap/test/pyout/aresult/ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc'


# time periods
cal_period = slice('2000-01-01', '2015-12-31')
val_period = slice('2016-01-01', '2016-12-31')
plot_period = slice('2016-09-03', '2030-10-13')

# tscale_sim dir
tscale_sim_dir = "/home/joel/sim/qmap/GR_data/sim/g3/"
grid = tscale_sim_dir.split('/')[-2]
grid = 'g3'
# ===============================================================================
# qmap_hor-plots.R settings
# ===============================================================================
# grid loop

indir = wd + '/' + grid

if not os.path.exists(indir):
    os.makedirs(indir)

lp = pd.read_csv(tscale_sim_dir + "/listpoints.txt")
mylon = lp.lon.mean()  # mean(lp.lon) # normally all lon are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre
mylat = lp.lat.mean()  # normally all lat are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre
tz = lp.tz.mean()
tscale_files = glob.glob(tscale_sim_dir + "/forcing/" + "meteoc*")

# run topoClim precprocessing to generate files corresponding to grid centre
outdir = indir
# just take first one for dissagregation as they are all scaled versions
# of one another - the signal should be fine but need to test. Incentive
# is to run this only once per grid = large efficiency gains
path_inp = tscale_files[0]
tclim.main(
    raw_dir,
    mylon,
    mylat,
    tz,
    nc_standard_clim,
    nc_standard_hist,
    cal_period,
    val_period,
    plot_period,
    path_inp,
    outdir)


# ===============================================================================
# code
# ===============================================================================


# loop over files
for tscale_file in range(len(tscale_files)):

    filename = tscale_files[tscale_file]
    sampID = filename.split('meteoc')[1].split('.csv')[0]

    indir2 = indir + '/s' + sampID
    indir_input = indir
    if not os.path.exists(indir2):
        os.makedirs(indir2)

    # station attributes
    slp = lp.slp[tscale_file]  # used in tclim_convert.py

    # linearly resample 3h era5 to 1h and cpture path as variable
    path_inp = resamp.main(tscale_files[tscale_file])

    root = indir
    grid = sampID

    # run qmap
    #cmd = ["Rscript", "qmap_hour_plots.R", indir2,indir_input,  path_inp]
    cmd = ["Rscript", "qmap_hour_plots.R", root, grid, path_inp]
    subprocess.check_output(cmd)

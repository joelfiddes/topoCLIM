import sys
import os
sys.path.append(os.getcwd()) 
sys.path.append('/home/caduff/src/topoCLIM/tclim3/lib/python3.7/site-packages')
import run_qmap_HPC
import pandas as pd
import glob
import logging
import resample_timeseries as resamp
import topoCLIM as tclim
from joblib import Parallel, delayed


num_cores = int(sys.argv[1]) # Must input number of cores

# make sure pacjaes are found
# https://stackoverflow.com/questions/39241032/how-to-import-a-local-python-module-when-using-the-sbatch-command-in-slurm


wd = "/home/joel/sim/qmap"
tscale_sim_dir = wd+ "/GR_data/sim/g4/"
indir = wd+"/topoclim_test_hpc/"

#===============================================================================
# hardcoded stuff
#===============================================================================
raw_dir= wd+'/raw_cordex/'


if not os.path.exists(indir):
	os.makedirs(indir)


# to clear logger: https://stackoverflow.com/questions/30861524/logging-basicconfig-not-creating-log-file-when-i-run-in-pycharm
for handler in logging.root.handlers[:]:
	logging.root.removeHandler(handler)

logging.basicConfig(level=logging.DEBUG, filename=indir+"/logfile",filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")



# standard calender examples (can be any file with standard cal) to interp non-standard cals to in hist and clim period
nc_standard_clim=raw_dir+'/aresult/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
nc_standard_hist=raw_dir+'/aresult/ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc' # exists?

# time periods
cal_period = slice('2000-01-01', '2015-12-31')
val_period = slice('2016-01-01', '2016-12-31')
plot_period = slice('2016-09-03', '2030-10-13')

# tscale_sim dir

grid = tscale_sim_dir.split('/')[-2]

#===============================================================================
# qmap_hor-plots.R settings
#===============================================================================
# grid loop

root = indir+'/'+ grid

if not os.path.exists(root):
	os.makedirs(root)

lp = pd.read_csv(tscale_sim_dir + "/listpoints.txt")
mylon = lp.lon.mean()#mean(lp.lon) # normally all lon are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre
mylat = lp.lat.mean()# normally all lat are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre
tz = lp.tz.mean()
tscale_files = glob.glob(tscale_sim_dir+"/forcing/"+ "meteoc*")

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/forcing/"+ "*1H.csv"):
	os.remove(f)

# rerun
tscale_files = glob.glob(tscale_sim_dir+"/forcing/"+ "meteoc*")

# run topoClim precprocessing to generate files corresponding to grid centre

path_inpt = tscale_files[0] # just take first one for dissagregation as they are all scaled versions of one another - the signal should be fine but need to test. Incentive is to run this only once per grid = large efficiency gains
path_inpt_1H = resamp.main(path_inpt)

logging.info("run topoclim")
tclim.main(raw_dir, mylon, mylat, tz, nc_standard_clim, nc_standard_hist, cal_period, val_period, plot_period, path_inpt_1H, root)

# delete path_inpt_1H sonot reprocessed
os.remove(path_inpt_1H)

# all jobs
tscale_files_sort = sorted(tscale_files)
print("running jobs: "+str(tscale_files_sort))

njobs=len(tscale_files_sort)

logging.info("starting " +str(njobs)+ " jobs")
Parallel(n_jobs=int(num_cores))(delayed(run_qmap_HPC.main)(root, tscale_file  ) for tscale_file in tscale_files_sort)

print("All cluster jobs complete!")

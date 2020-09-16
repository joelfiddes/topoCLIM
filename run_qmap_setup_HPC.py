import sys
import os
sys.path.append(os.getcwd()) 
import run_qmap_HPC
import pandas as pd
import glob

from joblib import Parallel, delayed

num_cores = int(sys.argv[1]) # Must input number of cores

# make sure pacjaes are found
# https://stackoverflow.com/questions/39241032/how-to-import-a-local-python-module-when-using-the-sbatch-command-in-slurm


#===============================================================================
# hardcoded stuff
#===============================================================================
wd = "/home/joel/sim/qmap"
raw_dir= wd+'/raw_cordex/'
indir = wd + "/topoclim_test_hpc/"
if not os.path.exists(indir):
	os.makedirs(indir)

# standard calender examples (can be any file with standard cal) to interp non-standard cals to in hist and clim period
nc_standard_clim=raw_dir+'/aresult/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
nc_standard_hist=raw_dir+'/aresult/ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc'

# time periods
cal_period = slice('2000-01-01', '2015-12-31')
val_period = slice('2016-01-01', '2016-12-31')
plot_period = slice('2016-09-03', '2030-10-13')

# tscale_sim dir
tscale_sim_dir = "/home/joel/sim/qmap/GR_data/sim/g3/"
grid = tscale_sim_dir.split('/')[-2]
grid='g3'
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

# run topoClim precprocessing to generate files corresponding to grid centre

#path_inp = tscale_files[0] # just take first one for dissagregation as they are all scaled versions of one another - the signal should be fine but need to test. Incentive is to run this only once per grid = large efficiency gains
#tclim.main(raw_dir, mylon, mylat, tz, nc_standard_clim, nc_standard_hist, cal_period, val_period, plot_period, path_inp, root)



# all jobs
simdirs = sorted(tscale_files)
print("running jobs: "+str(simdirs))

njobs=len(simdirs)


Parallel(n_jobs=int(num_cores))(delayed(run_qmap_HPC.main)(root, i+1, tscale_files[i]  ) for i in range(0,njobs))

print("All cluster jobs complete!")

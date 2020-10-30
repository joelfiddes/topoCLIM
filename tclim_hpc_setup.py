# a more evolved way of tclim' ing

"""
python3


Example:


Vars:


Details:


"""
import sys
import os
from joblib import Parallel, delayed
import logging
import tclim_src as tclim
import glob

wd=sys.argv[1] #'/home/joel/sim/qmap/topoclim_ch/'
raw_dir = sys.argv[2] # /home/caduff/sim/tclim/raw_cordex
tscale_sim_dir =sys.argv[3] 
num_cores=sys.argv[4] #10

#===============================================================================
# INPUT
#===============================================================================
nc_standard_clim=raw_dir+'/aresult/standard/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
nc_standard_hist=raw_dir+'/aresult/standard/ICHEC-EC-EARTH_historical_r12i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc'

# =========================================================================
#	Log / SETUP
# =========================================================================
if not os.path.exists(wd):
		os.makedirs(wd)

if not os.path.exists(wd+"/logs/"):
		os.makedirs(wd+"/logs/")

logfile = wd+ "/logs/logfile_setup"
if os.path.isfile(logfile) == True:
    os.remove(logfile)


# to clear logger: https://stackoverflow.com/questions/30861524/logging-basicconfig-not-creating-log-file-when-i-run-in-pycharm
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(level=logging.DEBUG, filename=logfile,filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")

logging.info("Run script = " + os.path.basename(__file__))


# =========================================================================
#	Kode
# =========================================================================
# list avaliable cordex
nc_complete = tclim.completeFiles(raw_dir) 
print(nc_complete)
logging.info(nc_complete)
# convert all cordex to standard calender - should be run once per domain (but is quick)
Parallel(n_jobs=int(num_cores))(delayed(tclim.calendarNinja)(nc,nc_standard_hist,nc_standard_clim) for nc in nc_complete)
print("CalendarNinja done!")

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1D.csv"):
	os.remove(f)

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1H.csv"):
	os.remove(f)


# find all era5 meteo files after cleanup
tscale_files = sorted(glob.glob(tscale_sim_dir+"/out/"+ "tscale*"))

logging.info("Number tscalefiles= "+str(len(tscale_files)))
logging.info("Number tscalefiles= "+str(len(tscale_files)))




















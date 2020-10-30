# a more evolved way of tclim' ing

"""
python3

This module takes post processed daily CORDEX downloads from esgf_post.py and 
produces point timeseries with standard calenders

Example:


Vars:


Details:


"""
import sys
import os
import glob
import xarray as xr
import calendar3 as cal3
import shutil
import subprocess
import pandas as pd
from tqdm import tqdm
import tclim_disagg
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import logging
import tclim_src as tclim

#grid= sys.argv[1]
wd=sys.argv[1]
tscale_sim_dir=sys.argv[2]
cordex_dir=sys.argv[3]
starti = sys.argv[4]
endi= sys.argv[5]
#===============================================================================
# INPUT
#===============================================================================
#wd = '/home/joel/sim/qmap/topoclim_ch/'
#tscale_sim_dir = "/home/joel/sim/qmap/ch_tmapp2/"
CORDEXPATH=cordex_dir+"/aresult/"
jobid = os.getenv('SLURM_ARRAY_TASK_ID')
# =========================================================================
#	Log
# =========================================================================
logfile = wd+ "/logs/logfile_qmap"+str(jobid)
if os.path.isfile(logfile) == True:
    os.remove(logfile)


# to clear logger: https://stackoverflow.com/questions/30861524/logging-basicconfig-not-creating-log-file-when-i-run-in-pycharm
for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)

logging.basicConfig(level=logging.DEBUG, filename=logfile,filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")
logging.info("Run script = " + os.path.basename(__file__))

#===============================================================================
# KODE
#===============================================================================

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1D.csv"):
	if os.path.exists(f):  
		os.remove(f)

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1H.csv"):
	if os.path.exists(f):  
		os.remove(f)

# get grid box
lp = pd.read_csv(tscale_sim_dir + "/listpoints.txt")

logging.info("Computing qmap files " + str(range(int(starti)-1,int(endi)) ) )

# find all era5 meteo files after cleanup
tscale_files = sorted(glob.glob(tscale_sim_dir+"/out/"+ "tscale*"))

# natural sorting https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

tscale_files.sort(key=natural_keys)

mytasks = range(int(starti)-1,int(endi))
for i in mytasks:
	tscale_file = tscale_files[i]
	logging.info("qmap " + tscale_file)
	print("qmap " + tscale_file)
	
	daily_obs = tclim.resamp_1D(tscale_file)
	sample =daily_obs.split('/')[-1].split('.')[0]
	cmd = ["Rscript", "qmap_hour_plots_daily.R", wd ,str(sample),  daily_obs, str(lp.lon[i]), str(lp.lat[i]), CORDEXPATH]
	subprocess.check_output(cmd)



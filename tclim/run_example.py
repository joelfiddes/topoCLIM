"""
python3

This module takes post processed daily CORDEX downloads from esgf_post.py and
produces point timeseries with standard calenders

Example:


Vars:


Details:


"""
import re
import sys
import os
import glob
import subprocess
import tclim_src as tclim
import tclim_disagg
import pandas as pd

wd = sys.argv[1]
starti = sys.argv[2]
endi = sys.argv[3]

# ===============================================================================
# INPUT
# ===============================================================================

# settings for example
# wd = "./"
# starti = 1
# endi = 1

tscale_sim_dir = wd+"/tscale/"
cordex_dir = wd+"/cordex"
CORDEXPATH = cordex_dir

# =========================================================================
#	Log
# =========================================================================
# logfile = wd+ "/logs/logfile_qmap"+str(jobid)
# if os.path.isfile(logfile) == True:
#     os.remove(logfile)


# # to clear logger: https://stackoverflow.com/questions/30861524/logging-basicconfig-not-creating-log-file-when-i-run-in-pycharm
# for handler in logging.root.handlers[:]:
#     logging.root.removeHandler(handler)

# logging.basicConfig(level=logging.DEBUG, filename=logfile,filemode="a+",format="%(asctime)-15s %(levelname)-8s %(message)s")
# logging.info("Run script = " + os.path.basename(__file__))

# ===============================================================================
# CODE
# ===============================================================================


# get grid box
lp = pd.read_csv(tscale_sim_dir + "/listpoints.txt")

print("Computing qmap files " + str(range(int(starti) - 1, int(endi))))

# find all tscale file excluding 1H and 1D ones
tscale_files = sorted(glob.glob(tscale_sim_dir + "/out/" + "tscale*"))
b = [item for item in tscale_files if '1H' not in item]
tscale_files = [item for item in b if '1D' not in item]

# natural sorting
# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [atoi(c) for c in re.split(r'(\d+)', text)]


tscale_files.sort(key=natural_keys)

mytasks = range(int(starti) - 1, int(endi))

for i in mytasks:

    tscale_file = tscale_files[i]

    print("qmap " + tscale_file)

    daily_obs = tclim.resamp_1D(tscale_file)
    sample = daily_obs.split('/')[-1].split('.')[0]

    if os.path.exists(wd + '/s' + sample + "/QSUCCESS"):
        print(tscale_file + " done!")
        continue

    print("qmap... ")
    cmd = ["Rscript", "../rsrc/qmap_hour_plots_daily_12.R", wd,
           str(sample), daily_obs, str(lp.lon[i]), str(lp.lat[i]), CORDEXPATH]
    subprocess.check_output(cmd)

    print("Aggregate results... ")
    cmd = ["Rscript", "../rsrc/aggregate_qmap_results.R", wd, str(sample)]
    subprocess.check_output(cmd)

    cmd = ["Rscript", "../rsrc/qmap_plots.R", wd, str(sample), daily_obs]
    subprocess.check_output(cmd)

    # cleanup
    tclim.findDelete(wd + "/s" + sample + "/aqmap_results", dir=True)

    # list all daily qmap files
    daily_cordex_files = glob.glob(wd + "/s" + sample + "/fsm/*Q.txt")
    # hourly obs
    hourly_obs = tclim.resamp_1H(tscale_file)
    # loop over with dissag routine

    print("Dissagregate time... ")
    for daily_cordex in daily_cordex_files:
        tclim_disagg.main(
            daily_cordex, hourly_obs, str(
                lp.lon[i]), str(
                lp.lat[i]), str(
                lp.tz[i]), str(
                    lp.slp[i]))

    meteofiles = (sorted(glob.glob(wd + "/s" + sample + "/fsm/*F.txt")))

    f = open(wd + '/s' + sample + "/QSUCCESS", "w")
    print("TopoCLIM finished!")


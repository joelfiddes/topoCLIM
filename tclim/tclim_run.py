
"""
    DESCRIPTION:
        This helper function downloads CORDEX data from the ESGF network. It requires 
        a user account to be set up first at:  
        https://esgf-data.dkrz.de/user/add/?next=http://esgf-data.dkrz.de/user/add/

    ARGS:
        wd (str): working directory of simulation, already exists
        tscale_dir (str): full path to location of TopoSCALE files (normally in wd)
        cordex_path (str): full path to downloaded and processed (esgf_post.py) CORDEX data

    RETURNS:
        NULL (files written to wd)
    


   """

import re
import sys
import os
import glob
import subprocess
import tclim_src as tsrc
import tclim_disagg as disagg
import pandas as pd

wd = sys.argv[1]
tscale_dir = sys.argv[2]
cordex_path = sys.argv[3]

# get grid box
lp = pd.read_csv(tscale_dir + "/listpoints.txt")

# find all tscale file excluding 1H and 1D ones
tscale_files = sorted(glob.glob(tscale_dir +  "/tscale*"))
b = [item for item in tscale_files if '1H' not in item]
tscale_files = [item for item in b if '1D' not in item]

# natural sorting
# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]


tscale_files.sort(key=natural_keys)
mytasks = range(0 , len(tscale_files) )

for i in mytasks:

    tscale_file = tscale_files[i]

    print("Tclim run: " + tscale_file)

    daily_obs = tsrc.resamp_1D(tscale_file)
    sample = daily_obs.split('/')[-1].split('.')[0]
    if os.path.exists(wd + '/s' + sample + "/QSUCCESS"):
        print(tscale_file + " done!")
        continue

    print("Quantile mapping... ")
    cmd = ["Rscript", "../rsrc/qmap_hour_plots_daily_12.R", wd,
           str(sample), daily_obs, str(lp.lon[i]), str(lp.lat[i]), cordex_path]
    subprocess.check_output(cmd)

    cmd = ["Rscript", "../rsrc/aggregate_qmap_results.R", wd, str(sample)]
    subprocess.check_output(cmd)

    # cleanup
    tsrc.findDelete(wd + "/s" + sample + "/aqmap_results", dir=True)

    # list all daily qmap files
    daily_cordex_files = glob.glob(wd + "/s" + sample + "/fsm/*Q.txt")
    # hourly obs
    hourly_obs = tsrc.resamp_1H(tscale_file)
    # loop over with dissag routine

    print("Dissagregate daily to hourly... ")
    for daily_cordex in daily_cordex_files:
        disagg.main(
            daily_cordex, hourly_obs, str(
                lp.lon[i]), str(
                lp.lat[i]), str(
                lp.tz[i]), str(
                    lp.slp[i]))


    f = open(wd + '/s' + sample + "/QSUCCESS", "w")

    print("TopoCLIM finished!")


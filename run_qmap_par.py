import pandas as pd
import glob
import resample_timeseries as resamp
import topoCLIM as tclim
import subprocess
import os
# asuume run:
# esgf_get.py
# esgf_post.py


def main(root, sample, path_inpt):
	#===============================================================================
	# qmap_hor-plots.R settings
	#===============================================================================
	# grid loop


	# tscale_files = glob.glob(tscale_sim_dir+"/forcing/"+ "meteoc*")
	# path_inp = tscale_files[0] # just take first one for dissagregation as they are all scaled versions of one another - the signal should be fine but need to test. Incentive is to run this only once per grid = large efficiency gains

	# # runs in setup once per grid
	# # tclim.main(mylon, mylat, tz, nc_standard_clim, nc_standard_hist, cal_period, val_period, plot_period, path_inp, outdir)


	# #===============================================================================
	# # code
	# #===============================================================================

	# # loop over files
	# tscale_file=tscale_files[file_index]
	# filename = tscale_files[tscale_file]
	# sampID = filename.split('meteoc')[1].split('.csv')[0]  


	# # linearly resample 3h era5 to 1h and cpture path as variable
	path_inpt_1H = resamp.main(path_inpt)
	# root = indir
	# sample= sampID

	# run qmap
	cmd = ["Rscript", "qmap_hour_plots.R", root ,str(sample),  path_inpt_1H]
	subprocess.check_output(cmd)


if __name__ == '__main__':
	root = sys.argv[1]
	sample = sys.argv[2]
	path_inpt =sys.argv[3]
	main(root, sample, path_inpt)

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


	# # linearly resample 3h era5 to 1h and cpture path as variable
	path_inpt_1H = resamp.main(path_inpt)

	# run qmap
	cmd = ["Rscript", "qmap_hour_plots.R", root ,str(sample),  path_inpt_1H]
	subprocess.check_output(cmd)


if __name__ == '__main__':
	root = sys.argv[1]
	sample = sys.argv[2]
	path_inpt =sys.argv[3]
	main(root, sample, path_inpt)

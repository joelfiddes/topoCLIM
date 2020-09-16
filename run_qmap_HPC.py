
import subprocess
import os

# asuume run:
# esgf_get.py
# esgf_post.py


def main(root, path_inpt_1H):

	# wd = "/home/joel/sim/qmap"
	# raw_dir= wd+'/raw_cordex/'
	# indir = wd + "/topoclim_test_hpc/"
	# if not os.path.exists(indir):
	# 	os.makedirs(indir)

	# # standard calender examples (can be any file with standard cal) to interp non-standard cals to in hist and clim period
	# nc_standard_clim=raw_dir+'/aresult/ICHEC-EC-EARTH_rcp85_r12i1p1_CLMcom-CCLM5-0-6_v1__TS.nc_TS_ALL_ll.nc'
	# nc_standard_hist=raw_dir+'/aresult/ICHEC-EC-EARTH_historical_r1i1p1_KNMI-RACMO22E_v1__TS.nc_TS_ALL_ll.nc'

	# # time periods
	# cal_period = slice('2000-01-01', '2015-12-31')
	# val_period = slice('2016-01-01', '2016-12-31')
	# plot_period = slice('2016-09-03', '2030-10-13')

	# # tscale_sim dir
	# tscale_sim_dir = "/home/joel/sim/qmap/GR_data/sim/g3/"
	# grid = tscale_sim_dir.split('/')[-2]
	# grid='g3'

	# lp = pd.read_csv(tscale_sim_dir + "/listpoints.txt")
	# mylon = lp.lon.mean()#mean(lp.lon) # normally all lon are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre
	# mylat = lp.lat.mean()# normally all lat are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre
	# tz = lp.tz.mean()

	# print(sample)
	# # # linearly resample 3h era5 to 1h and cpture path as variable
	# path_inpt_1H = resamp.main(path_inpt)

	# # just take first one for dissagregation as they are all scaled versions of one another - the signal should be fine but need to test. Incentive is to run this only once per grid = large efficiency gains
	# tclim.main(raw_dir, mylon, mylat, tz, nc_standard_clim, nc_standard_hist, cal_period, val_period, plot_period, path_inpt_1H, root)

	# run qmap

	sample =path_inpt_1H.split('/')[-1].split('.')[0]
	logging.info("Running qmap job "+ sample )
	cmd = ["Rscript", "qmap_hour_plots.R", root ,str(sample),  path_inpt_1H]
	subprocess.check_output(cmd)
	logging.info("Qmap job "+ sample + "complete") 

if __name__ == '__main__':
	root = sys.argv[1]
	path_inpt =sys.argv[2]
	main(root, path_inpt)

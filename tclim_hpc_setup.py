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

wd=sys.argv[1] #'/home/joel/sim/qmap/topoclim_ch/'
raw_dir = sys.argv[2] # /home/caduff/sim/tclim/raw_cordex
num_cores=sys.argv[3] #10

#===============================================================================
# INPUT
#===============================================================================
raw_dir= wd+'/raw_cordex/'
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

# convert all cordex to standard calender - should be run once per domain (but is quick)
Parallel(n_jobs=int(num_cores))(delayed(tclim.calendarNinja)(nc,nc_standard_hist,nc_standard_clim) for nc in nc_complete)





























# find all era5 meteo files
tscale_files =sorted(glob.glob(tscale_sim_dir+"/out/"+ "tscale*")  )

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1D.csv"):
	os.remove(f)

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1H.csv"):
	os.remove(f)

# get grid box
lp = pd.read_csv(tscale_sim_dir + "/listpoints.txt")
# this doesnt seem to be era5 grid centres - check
lon = lp.lon#mean(lp.lon) # normally all lon are the same (grid centre), however recent version tsub allows the position to be weight by pixel positions, mean() then gets back to grid centre, BUT needs weighting by sample memberes!
lat = lp.lat
tz = lp.tz

# rerun after cleanup
tscale_files = sorted(glob.glob(tscale_sim_dir+"/out/"+ "tscale*"))
Parallel(n_jobs=int(num_cores))(delayed(tclim_main)(tscale_files[i], lon[i], lat[i], tz[i], lp.slp[i]) for i in tqdm(range( len(tscale_files))) )

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1D.csv"):
	os.remove(f)

# clean up old resamples
for f in glob.glob(tscale_sim_dir+"/out/"+ "*1H.csv"):
	os.remove(f)

fsmfiles = (sorted(glob.glob(root+"/*/output/*.txt")))
plot_hs(fsmfiles)


tsub_root = "/home/joel/sim/qmap/ch_tmapp2"
root = "/home/joel/sim/qmap/topoclim_ch"
nsims=2100 # 2100
col=2
nclust = 100
# compile timeseries means HIST


start='1980-01-01' 
end = '2005-12-31'
scenario="HIST" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,"/home/joel/src/topoCLIM/spatialize.R" ,tsub_root ,meanVar, str(nclust), outname]
subprocess.check_output(cmd)

# compile timeseries meand RCP26 NEAR
start='2030-01-01' 
end = '2050-12-31'
scenario="RCP26" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,"/home/joel/src/topoCLIM/spatialize.R" ,tsub_root ,meanVar, str(nclust), outname]
subprocess.check_output(cmd)


# compile timeseries meand RCP26 FAR
start='2080-01-01' 
end = '2099-12-31'
scenario="RCP26" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,"/home/joel/src/topoCLIM/spatialize.R" ,tsub_root ,meanVar, str(nclust), outname]
subprocess.check_output(cmd)


# compile timeseries meand RCP85 NEAR
root = "/home/joel/sim/qmap/topoclim_ch"

start='2030-01-01' 
end = '2050-12-31'
scenario="RCP85" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,"/home/joel/src/topoCLIM/spatialize.R" ,tsub_root ,meanVar, str(nclust), outname]
subprocess.check_output(cmd)


# compile timeseries meand RCP85 FAR
start='2080-01-01' 
end = '2099-12-31'
scenario="RCP85" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,"/home/joel/src/topoCLIM/spatialize.R" ,tsub_root ,meanVar, str(nclust), outname]
subprocess.check_output(cmd)

# converts standard output to FSm or...
#qfiles = tqdm(sorted(glob.glob(root+"/*/fsm/*Q_H.txt")))
#qfiles = tqdm(sorted(glob.glob(root+"/smeteoc6_1D/fsm/*Q_H.txt")))
#for qfile in tqdm(qfiles):
#Parallel(n_jobs=int(num_cores))(delayed(met2fsm_parallel)(qfile,lp) for qfile in qfiles)
# cleanup _HOURLY.txt
#findDelete(root+"/*/fsm/*Q_H.txt")

# Simulate FSm
#meteofiles = tqdm(sorted(glob.glob(root+"/*/fsm/meteo/*F.txt")))
#meteofiles = tqdm(sorted(glob.glob(root+"/smeteoc6_1D/fsm/meteo/*F.txt")))
#os.chdir(indir)

# can 
#Parallel(n_jobs=int(1))(delayed(fsm_sim)(meteofile,namelist,fsmexepath) for meteofile in meteofiles)

# At present a simple fixed output format is used. The output text file has 10 columns:

# | Variable | Units  | Description       |
# |----------|--------|-------------------|
# | year     | years  | Year              |
# | month    | months | Month of the year |
# | day      | days   | Day of the month  |
# | hour     | hours  | Hour of the day   |
# | alb      | -      | Effective albedo  |
# | Rof      | kg m<sup>-2</sup> | Cumulated runoff from snow    |
# | snd      | m      | Average snow depth                       |
# | SWE      | kg m<sup>-2</sup> | Average snow water equivalent |
# | Tsf      | &deg;C | Average surface temperature              |
# | Tsl      | &deg;C | Average soil temperature at 20 cm depth  |






# map

# max filter = 4 (HS)
# max filter = 


# zmax=400
# maxfilter=2000
# var=3
# outname="swe"
# myroot='/home/joel/sim/qmap/topoclim_test_hpc'
# for grid in range(6):
# 	print(grid+1)

# 	root = '/home/joel/sim/qmap/topoclim_test_hpc/g'+ str(grid +1)
# 	tscale_sim_dir = '/home/joel/sim/qmap/GR_data/sim/g' +str(grid +1)
# 	spatialfsm(root,var, maxfilter)
# 	plot_map(root, tscale_sim_dir+"/landform.tif", str(zmax))	
	
# compile_map(myroot, outname, str(zmax))



# zmax=400
# maxfilter=2000
# var=3
# outname="swe"
# myroot='/home/joel/sim/qmap/GR_data'

# for grid in range(6):
# 	tscale_sim_dir = '/home/joel/sim/qmap/GR_data/sim/g' +str(grid +1)
# 	fsm_path='/home/joel/sim/qmap/GR_data/sim/g'+str(grid +1)+'/out/FSM'
# 	spatialfsm_era5(fsm_path,var)
# 	plot_map_era5(fsm_path, tscale_sim_dir+"/landform.tif", str(zmax))

# compile_map_era5(myroot, outname, str(zmax))
















# x5 compression 50GB -> 10GB
#findCompress(root + "/**/meteo")  
#findDelete(root+"/**/meteo", dir=True)

# def plot_fsm(file, col):
# 	df =pd.read_csv(file, delim_whitespace=True, parse_dates=[[0,1,2]], header=None) 
# 	df.set_index(df.iloc[:,0], inplace=True)  
# 	df.drop(df.columns[[0]], axis=1, inplace=True )  
# 	swe=df.iloc[:,col]
# 	plot(swe)


# def experiments():

# 	file="/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/out.txt"  
# 	file="/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/CNRM-CERFACS-CNRM-CM5_SMHI-RCA4_HIST_F.txt"

# 	file="/home/joel/sim/qmap/GR_data/sim/g3/out/FSM/fsm001.txt_00.txt"


# 	file="/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/out.txt" 


# 	fsmfiles = (sorted(glob.glob(root+"/*/fsm/output/*11111.txt")))

# 	root="//home/joel/sim/qmap/GR_data/sim/g3/out/FSM_config31/"
# 	fsmfiles = (sorted(glob.glob(root+"/fsm*")))






# 	file = '/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM5-0-6_HIST_Q.txt' 
# 	q1 =pd.read_csv(file, parse_dates=True)

# 	file="/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM5-0-6_HIST_Q_H.txt"
# 	hr1 =pd.read_csv(file, parse_dates=True)
# 	file="/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM5-0-6_HIST_Q_HOURLY.txt"
# 	hr2 =pd.read_csv(file, parse_dates=True)


# 	file= "/home/joel/sim/qmap/GR_data/forcing/fsm006.txt"
# 	era5 =pd.read_csv(file, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)

# 	# to tes
# 	filestotest = glob.glob("/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/meteo/*_F.txt")

# 	for f in filestotest:
# 		clim =pd.read_csv(f,delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
# 		print(clim.iloc[:,6].mean()  )


# 	plot(era5.iloc[:,2], clim.iloc[:,2])

# 	era5.iloc[:,6].plot()
# 	clim.iloc[:,2].plot()
# 	plt.show()

# 	plt.scatter(era5.iloc[:,6] , clim.iloc[:,2])


# 	# configs are good 00 and 31 more or less the same

# 	# test all forcings
# 	col=9
# 	file="/home/joel/sim/qmap/GR_data/forcing/fsm006.txt"
# 	df =pd.read_csv(file, delim_whitespace=True, parse_dates=[[0,1,2,3]], header=None)
# 	df.set_index(df.iloc[:,0], inplace=True)  
# 	df.drop(df.columns[[0]], axis=1, inplace=True )
# 	era51h =df.resample('1H').interpolate()
# 	print(era5.iloc[:,col].mean()  )


# 	filestotest = glob.glob("/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/meteo/*HIST_F.txt")

# 	climmean=[]
# 	for f in filestotest:
# 		clim =pd.read_csv(f,delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
# 		climmean.append(clim.iloc[:,col].mean() )
# 	print(np.mean(climmean  ))

# 	col=6
# 	filestotest = (sorted(glob.glob(root+"/smeteoc6_1D/fsm/output/*11111.txt")))
# 	for f in filestotest:
# 		df =pd.read_csv(f,delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
# 		df.set_index(df.iloc[:,0], inplace=True)  
# 		df.drop(df.columns[[0]], axis=1, inplace=True )  
# 		swe=df.iloc[:,2]
# 		swe.plot(title="ID:"+ f  )
# 		plt.show()


# 	# tests show sw is half and ta has -3deg bias, by simple bias correction results are much better - what is wrong with qmap /diaagg?
# 	# now test how daily data looks	


# 	file = '/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM5-0-6_HIST_Q.txt' 
# 	q1 =pd.read_csv(file, parse_dates=True)
# 	file="/home/joel/sim/qmap/topoclim_test_hpc/g3/smeteoc6_1D/fsm/CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM5-0-6_HIST_Q_H.txt"
# 	hr1 =pd.read_csv(file, parse_dates=True)
# 	file= "/home/joel/sim/qmap/GR_data/forcing/fsm006.txt"
# 	era5 =pd.read_csv(file, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)


# In [78]: q1.ISWR.mean()                                                                                                                                                             
# Out[78]: 152.95145001013992                                                                                                                                                         
# In [79]: hr1.ISWR.mean()                                                                                                                                                            
# Out[79]: 49.450357362678844                                                                                                                                                         
# In [82]: era5.iloc[:,2].mean()                                                                                                                                                      
# Out[82]: 119.06607230607231    

# In [83]: hr1.TA.mean()                                                                                                                                                              
# Out[83]: 271.824690732103                                                                                                                                                           
# In [84]: q1.TA.mean()                                                                                                                                                               
# Out[84]: 271.748326911377 
# In [87]: era5.iloc[:,6].mean()                                                                                                                                                      
# Out[87]: 274.4854905229906     




#EVALUATION

# grid = 'g4'
# sample='5'
# file = "/home/joel/sim/qmap/topoclim_test_hpc/"+grid+"/smeteoc"+sample+"_1D/fsm/ICHEC-EC-EARTH_CLMcom-CCLM5-0-6_HIST_Q.txt" 
# q1 =pd.read_csv(file, parse_dates=True)
# q1 =pd.read_csv(file, parse_dates=True)
# q1.set_index(q1.iloc[:,0], inplace=True)
# q1.index = pd.to_datetime(q1.index)


# file="/home/joel/sim/qmap/topoclim_test_hpc/"+grid+"/smeteoc"+sample+"_1D/fsm/CNRM-CERFACS-CNRM-CM5_CLMcom-CCLM5-0-6_HIST_Q_H.txt"

# hr1 =pd.read_csv(file, parse_dates=True)
# hr1.set_index(hr1.iloc[:,0], inplace=True)
# hr1.index = pd.to_datetime(hr1.index)

# file= "/home/joel/sim/qmap/GR_data/sim/g4/forcing/meteoc"+sample+"_1H.csv"
# era5 =pd.read_csv(file, parse_dates=True)
# era5.set_index(era5.iloc[:,0], inplace=True)
# era5.index = pd.to_datetime(era5.index)

# q1.ISWR.mean() 
# hr1.ISWR.resample('1d').mean().mean()
# era5.ISWR.mean() 


# print(q1.TA.mean())
# print(hr1.TA.resample('1d').mean().mean())  
# print(era5.TA.mean()  )


# In [78]: q1.ISWR.mean()                                                                                                                                                             
# Out[78]: 152.95145001013992                                                                                                                                                         
# In [79]: hr1.ISWR.resample('1d').mean()                                                                                                                                                            
# Out[79]: 49.450357362678844                                                                                                                                                         
# In [82]: era5.iloc[:,2].mean()                                                                                                                                                      
# Out[82]: 119.06607230607231    

# In [83]: hr1.TA.mean()                                                                                                                                                              
# Out[83]: 271.824690732103                                                                                                                                                           
# In [84]: q1.TA.mean()                                                                                                                                                               
# Out[84]: 271.748326911377 
# In [87]: era5.iloc[:,6].mean()                                                                                                                                                      
# Out[87]: 274.4854905229906     
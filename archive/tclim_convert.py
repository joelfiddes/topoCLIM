# converts standard output to FSm or...
import glob
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
slope=0
fsm_out=True
grid="g4"
#indir="/home/joel/sim/qmap/topoclim/fsm/"
indir="/home/joel/sim/qmap/topoclim_test_hpc"
os.chdir(indir)   
#outdir=indir+"/meteo/"
qfiles = glob.glob(indir + "/"+grid+"/*/fsm/*QMAP.txt")




for qfile in qfiles:
	print(qfile)
	qdat= pd.read_csv(qfile, index_col=0, parse_dates=True)
	outdir = os.path.split(qfile)[0]  + "/meteo/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)


	#======================================================================================
	# # partition prate to rain snow (mm/hr)
	#======================================================================================
			
	lowthresh=272.15
	highthresh = 274.15
	d = {'prate': qdat.PINT, 'ta': qdat.TA}
	df = pd.DataFrame(data=d)
	snow = df.prate.where(df.ta<lowthresh) 
	rain=df.prate.where(df.ta>=highthresh) 

	mix1S = df.prate.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mix1T = df.ta.where((df.ta >= lowthresh) & (df.ta<=highthresh), inplace = False)
	mixSno=(highthresh - mix1T) / (highthresh-lowthresh)
	mixRain=1-mixSno
	addSnow=mix1S*mixSno
	addRain=mix1S*mixRain

	# nas to 0
	snow[np.isnan(snow)] = 0 	
	rain[np.isnan(rain)] = 0 
	addRain[np.isnan(addRain)] = 0 
	addSnow[np.isnan(addSnow)] = 0 

	#======================================================================================
	# # linearly reduce snow to zero in steep slopes
	#if steepSnowReduce=="TRUE": # make this an option if need that in future
	#======================================================================================

	snowSMIN=30.
	snowSMAX=80.
	slope=slope

	k= (snowSMAX-slope)/(snowSMAX - snowSMIN)

	if slope<snowSMIN:
		k=1
	if slope>snowSMAX:
		k=0

	snowTot=(snow+addSnow) * k
	rainTot=rain + addRain



	if fsm_out is True:
		outname1=qfile.split('/')[-1]  
		outname=outname1.split('__TS.nc')[0]


		#Filter TA for absolute 0 vals
		
		qdat.TA[qdat.TA<220]=np.nan
		qdat.TA = qdat.TA.ffill()
		qdat.P[qdat.P<10000]=np.nan
		qdat.P = qdat.P.ffill()
		dates=qdat.index
		df_fsm= pd.DataFrame({	
		 				"year": dates.year, 
		 				"month": dates.month, 
						"day": dates.day, 
						"hour": dates.hour,
						"ISWR":qdat.ISWR, 
						"ILWR":qdat.ILWR, 
						"Sf":snowTot/(60.*60.), # prate in mm/hr to kgm2/s
						"Rf":rainTot/(60.*60.), # prate in mm/hr to kgm2/s
						"TA":qdat.TA, 
						"RH":qdat.RH,#*0.01, #meteoio 0-1
						"VW":qdat.VW,
						"P":qdat.P,
						
						
						})

		df_fsm.to_csv(path_or_buf=outdir+outname+"_FSM.txt" ,na_rep=-999,float_format='%.8f', header=False, sep='\t', index=False, 
			columns=['year','month','day', 'hour', 'ISWR', 'ILWR', 'Sf', 'Rf', 'TA', 'RH', 'VW', 'P'])


#========================================================================
# Sim - make sure cd into dir
#========================================================================

import numpy as np
import os
import sys


namelist="/home/joel/sim/qmap/topoclim/fsm/nlst_qmap.txt"
#os.system('./compil.sh')

#METEOFILES  = glob.glob(outdir + "*FSM.txt")
METEOFILES = glob.glob(indir + "/"+grid+"/*/fsm/meteo/*FSM.txt")




os.chdir(indir)
for METEOFILE in tqdm(METEOFILES):
	print(METEOFILE)

	METEOFILENAME = os.path.split(METEOFILE)[1]
	METEOFILEPATH = os.path.split(METEOFILE)[0]
	FSMPATH = os.path.split(METEOFILEPATH)[0] 
	try:
		os.mkdir(FSMPATH+'/output')
	except:
    	pass


	#for n in 31: #range(32):
	n=31

    config = np.binary_repr(n, width=5)
    print('Running FSM configuration ',config,n)
    f = open('nlst.txt', 'w')
    out_file = 'out.txt'
    with open(namelist) as file:
        for line in file:
            f.write(line)
            if 'config' in line:
                f.write('  nconfig = '+str(n)+'\n')
            if 'drive' in line:
                f.write('  met_file = ' +"'./meteo/"+METEOFILENAME+"'"+'\n')
            if 'out_file' in line:
                out_file = line.rsplit()[-1]
            out_name = out_file.replace('.txt','')
    f.close()

    os.system("cp FSM " +FSMPATH)
    os.system("cp nlst.txt " +FSMPATH)
    os.chdir(FSMPATH)
    os.system('./FSM < nlst.txt')
    save_file = FSMPATH + '/output/'+METEOFILENAME+ out_name+'_'+config+'.txt'
    os.system('mv '+out_file+' '+save_file)
		#os.system('rm nlst.txt')


#========================================================================
# Plot models
#========================================================================


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import matplotlib

var=2 # HS
#var=3 # swe

# swe
def plot(obs):
    plt.figure()
    ax = plt.gca()
    obs.plot()
    plt.legend()
    plt.show()



stats_hist=pd.DataFrame()
stats_rcp=pd.DataFrame()

ncfiles = glob.glob("/home/joel/sim/qmap/topoclim/fsm/output/"+ "*HIST*")
swe_hist=pd.DataFrame()

for nc in ncfiles:
	name1 = nc.split('/')[-1]
	name = name1.split('.')[0]
	df_obs= pd.read_csv(nc, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
	df_obs.set_index(df_obs.iloc[:,0], inplace=True)  
	df_obs.drop(df_obs.columns[[0]], axis=1, inplace=True )      
	swe_hist[name]=df_obs.iloc[:,var]
	
#plot(swe_hist.mean(axis=1 )) 
stats_hist['histmean'] =swe_hist.mean(axis=1 )
stats_hist['histstd'] =swe_hist.std(axis=1 )

ncfiles = glob.glob("/home/joel/sim/qmap/topoclim/fsm/output/"+ "*RCP26*")
swe_rcp26=pd.DataFrame()
for nc in ncfiles:
	name1 = nc.split('/')[-1]
	name = name1.split('.')[0]
	df_obs= pd.read_csv(nc, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
	df_obs.set_index(df_obs.iloc[:,0], inplace=True)  
	df_obs.drop(df_obs.columns[[0]], axis=1, inplace=True )  
	swe_rcp26[name]=df_obs.iloc[:,var]

#plot(swe_rcp26.mean(axis=1 )) 
stats_rcp['rcp26mean'] =swe_rcp26.mean(axis=1 )
stats_rcp['rcp26std'] =swe_rcp26.std(axis=1 )

ncfiles = glob.glob("/home/joel/sim/qmap/topoclim/fsm/output/"+ "*RCP85*")
swe_rcp85=pd.DataFrame()
for nc in ncfiles:
	name1 = nc.split('/')[-1]
	name = name1.split('.')[0]
	df_obs= pd.read_csv(nc, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
	df_obs.set_index(df_obs.iloc[:,0], inplace=True)  
	df_obs.drop(df_obs.columns[[0]], axis=1, inplace=True )  
	swe_rcp85[name]=df_obs.iloc[:,var]

#plot(swe_rcp85.mean(axis=1 )) 
stats_rcp['rcp85mean'] =swe_rcp85.mean(axis=1 )
stats_rcp['rcp85std'] =swe_rcp85.std(axis=1 )




hist_mean_year = stats_hist.histmean.resample('A').mean()
rcp26_mean_year = stats_rcp.rcp26mean.resample('A').mean()
rcp85_mean_year = stats_rcp.rcp85mean.resample('A').mean()

hist_std_year = stats_hist.histstd.resample('A').mean()
rcp26_std_year = stats_rcp.rcp26std.resample('A').mean()
rcp85_std_year = stats_rcp.rcp85std.resample('A').mean()



#========================================================================
# Sim ERA5 obs
#========================================================================

obsfile = "/home/joel/sim/qmap/wfj_long_1H.csv"
odat= pd.read_csv(obsfile, index_col=0, parse_dates=True)

dates=odat.index
df_fsm= pd.DataFrame({	
		 				"year": dates.year, 
		 				"month": dates.month, 
						"day": dates.day, 
						"hour": dates.hour,
						"ISWR":odat.ISWR, 
						"ILWR":odat.ILWR, 
						"Sf":odat.Sf/(60.*60.), # prate in mm/hr to kgm2/s
						"Rf":odat.Rf/(60.*60.), # prate in mm/hr to kgm2/s
						"TA":odat.TA, 
						"RH":odat.RH*100,#*0.01, #meteoio 0-1
						"VW":odat.VW,
						"P":odat.P,
						
						
						})

df_fsm.to_csv(path_or_buf=outdir+"/OBS_FSM.txt" ,na_rep=-999,float_format='%.8f', header=False, sep='\t', index=False, 
			columns=['year','month','day', 'hour', 'ISWR', 'ILWR', 'Sf', 'Rf', 'TA', 'RH', 'VW', 'P'])




import numpy as np
import os
import sys


namelist="/home/joel/sim/qmap/topoclim_test/fsm/nlst_qmap.txt"
#os.system('./compil.sh')

METEOFILES  = glob.glob(outdir + "OBS_FSM.txt")



try:
    os.mkdir('output')
except:
    pass

for METEOFILE in METEOFILES:
	print(METEOFILE)
	#for n in 31: #range(32):
	n=31
    config = np.binary_repr(n, width=5)
    print('Running FSM configuration ',config,n)
    f = open('nlst.txt', 'w')
    out_file = 'out.txt'
    with open(namelist) as file:
        for line in file:
            f.write(line)
            if 'config' in line:
                f.write('  nconfig = '+str(n)+'\n')
            if 'drive' in line:
            	METEOFILENAME=METEOFILE.split(indir+'/meteo/')[1]
                f.write('  met_file = ' +"'./meteo/"+METEOFILENAME+"'"+'\n')
            if 'out_file' in line:
                out_file = line.rsplit()[-1]
            out_name = out_file.replace('.txt','')
    f.close()
    os.system('./FSM < nlst.txt')
    save_file = 'output/'+METEOFILENAME+ out_name+'_'+config+'.txt'
    os.system('mv '+out_file+' '+save_file)
		#os.system('rm nlst.txt')

#plot gcos obs - no good as cant get mean annual values
gcos_dat="//home/joel/sim/qmap/wfj_GCOS.csv"
gcosdat= pd.read_csv(gcos_dat, index_col=1, parse_dates=True)



# get header line from smet
import re
pattern = "fields"

file = open("/home/joel/data/wfj_optimal/WFJ_optimaldataset_v8.smet", "r")
for line in file:
    if re.search(pattern, line):
        hdr = line

# silly smet hdr parsing
myhdr = hdr.split('=')[1] 
hdrout= myhdr[1:-1]



# read file
wfj_optim= pd.read_csv("/home/joel/data/wfj_optimal/WFJ_optimaldataset_v8.smet", index_col=0, parse_dates=True, skiprows=21, sep=r'\s{2,}' ,header=None, na_values=-999)
wfj_optim.columns= hdrout.split(' ')[1:]   

# aggregate to year
optim_year = wfj_optim.resample('A').mean()    

# swe
def plot(obs):
    plt.figure()
    ax = plt.gca()
    obs.plot()
    plt.legend()
    plt.show()

ncfiles = glob.glob("/home/joel/sim/qmap/topoclim_test/fsm/output/"+ "*OBS*")
for nc in ncfiles:

	df_obs= pd.read_csv(nc, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
	df_obs.set_index(df_obs.iloc[:,0], inplace=True)  
	df_obs.drop(df_obs.columns[[0]], axis=1, inplace=True )      
	era5=df_obs.iloc[:,var]

# annual values
era5_year = era5.resample('A').mean()
era5_year.name='era5'


# obs val

#========================================================================
# Plot
#========================================================================
plt.figure()
ax = plt.gca()

#era5_year.plot(color='grey') something wri=ong with this as biased where qmapped datra that depends on this is good!
optim_year['HS'][optim_year['HS']<0] = np.nan
wfj_HS = optim_year['HS'][0:-1] # drop incomplete 2019
wfj_HS.name="HS_obs" 
wfj_HS.plot(color='red')

plt.fill_between(hist_mean_year.index, hist_mean_year-hist_std_year , hist_mean_year+hist_std_year , alpha=0.2)
hist_mean_year.plot()

plt.fill_between(rcp85_mean_year.index, rcp85_mean_year-rcp85_std_year, rcp85_mean_year+rcp85_std_year, alpha=0.2)
rcp85_mean_year.plot()

plt.fill_between(rcp26_mean_year.index, rcp26_mean_year-rcp26_std_year, rcp26_mean_year+rcp26_std_year, alpha=0.2)
rcp26_mean_year.plot()
wfj_HS.plot(color='red')
plt.axhline(hist_mean_year.mean(), color='grey')
ax.set_xlabel('year')
ax.set_ylabel('mean annual snow height (m)')

plt.legend()
plt.show()















#HS
def plot(obs):
    plt.figure()
    ax = plt.gca()
    obs.plot()
    plt.legend()
    plt.show()

ncfiles = glob.glob('/home/joel/sim/qmap/topoclim/fsm/output/'+ "*.txt")
for nc in ncfiles:

	df_obs= pd.read_csv(nc, index_col=(0,1,2), delim_whitespace=True)
	plot(df_obs.iloc[:,2])

# TSS
def plot(obs):
    plt.figure()
    ax = plt.gca()
    obs.plot()
    plt.legend()
    plt.show()

ncfiles = glob.glob('/home/joel/sim/qmap/fsm/output1/'+ "*.txt")
for nc in ncfiles:

	df_obs= pd.read_csv(nc, index_col=(0,1,2), delim_whitespace=True)
	plot(df_obs.iloc[:,4])




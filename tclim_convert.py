# converts standard output to FSm or...
import glob
import pandas as pd
import numpy as np
slope=0
fsm_out=True
indir="/home/joel/sim/qmap/topoclim/"
outdir=indir+"/aqmap_results/"
qfiles = glob.glob(outdir + "*QMAP.txt")

for qfile in qfiles:
	qdat= pd.read_csv(qfile, index_col=0, parse_dates=True)


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
		outname=qfile+"_FSM.txt"
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

		df_fsm.to_csv(path_or_buf=outname ,na_rep=-999,float_format='%.8f', header=False, sep='\t', index=False, 
			columns=['year','month','day', 'hour', 'ISWR', 'ILWR', 'Sf', 'Rf', 'TA', 'RH', 'VW', 'P'])


#========================================================================
# Sim
#========================================================================

import numpy as np
import os
import sys

namelist = sys.argv[1]
#os.system('./compil.sh')



try:
    os.mkdir('output')
except:
    pass

for n in range(32):
    config = np.binary_repr(n, width=5)
    print('Running FSM configuration ',config,n)
    f = open('nlst.txt', 'w')
    out_file = 'out.txt'
    with open(namelist) as file:
        for line in file:
            f.write(line)
            if 'config' in line:
                f.write('  nconfig = '+str(n)+'\n')
            if 'out_file' in line:
                out_file = line.rsplit()[-1]
            out_name = out_file.replace('.txt','')
    f.close()
    os.system('./FSM < nlst.txt')
    save_file = 'output/'+out_name+'_'+config+'.txt'
    os.system('mv '+out_file+' '+save_file)
#os.system('rm nlst.txt')


#========================================================================
# Plot
#========================================================================




import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob
import matplotlib



# swe
def plot(obs):
    plt.figure()
    ax = plt.gca()
    obs.plot()
    plt.legend()
    plt.show()

ncfiles = glob.glob(outdir+ "*FSM.txt")
for nc in ncfiles:

	df_obs= pd.read_csv(nc, index_col=(0,1,2), delim_whitespace=True)
	plot(df_obs.iloc[:,3])




#HS
def plot(obs):
    plt.figure()
    ax = plt.gca()
    obs.plot()
    plt.legend()
    plt.show()

ncfiles = glob.glob('/home/joel/sim/qmap/fsm/output1/'+ "*.txt")
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
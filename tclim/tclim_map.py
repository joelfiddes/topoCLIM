import sys
import os
import glob
import subprocess
import pandas as pd
from tqdm import tqdm
import numpy as np

tsub_root = "/home/caduff/sim//ch_tmapp_50"
root = "/home/caduff/sim/tclim_ch"
rcode = "/home/caduff/src/topoCLIM/spatialize.R"
nsims=1050 # 2100
col=2
nclust = 50

tsub_root = "/home/caduff/sim//ch_tmapp_100"
root = "/home/caduff/sim/tclim_ch2"
rcode = "/home/caduff/src/tmapp2/rsrc/spatialize.R"
nsims=2100 # 2100
col=2
nclust = 100
Ngrid=21

def timeseries_means(root, nsims, col, start, end, scenario):
	print (scenario)
	print (start)
	import pandas as pd
	mean_ts=[]
	sd_ts=[]
	for ID in (range(nsims)):
		print(ID)
		filenames=root + "/stscale_"+str(ID+1)+"_1D"+"/output/*_"+scenario+"_Q_F.txt.txt"
		files = glob.glob(filenames)
		models=[]
		models_sd=[]
		for f in (files): # need to handle multiple models here
			try:
				df =pd.read_csv(f, delim_whitespace=True, parse_dates=[[0,1,2]], header=None)
			except:
				print(f + " failed read")
				
			df.set_index(df.iloc[:,0], inplace=True)  
			df.drop(df.columns[[0]], axis=1, inplace=True )  
			swe=df.iloc[:,col]
			try:
				swemean = swe[slice(start,end)].mean()
				swesd = swe[slice(start,end)].std()
				#mask = (df['date'].dt.month == 6) 
				#df.loc[mask]
				#df = pd.DataFrame(index=pd.date_range(dt.datetime(2015,1,1), dt.datetime(2015,5,1)))
				#df.loc[df.index.month==8] get august of every year
			except:
				print(f)
				print(str(ID)+ " failed")
			models.append(swemean)
			models_sd.append(swesd)
		mean_ts.append(np.mean(models)	)
		sd_ts.append(np.mean(models_sd)	)
	pd.Series(mean_ts).to_csv(root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv",header=False, float_format='%.3f') 
	pd.Series(sd_ts).to_csv(root+"/sd_ts_"+str(col)+"_"+scenario+start+end+".csv",header=False, float_format='%.3f') 







# compile timeseries means HIST


start='1980-01-01' 
end = '2005-12-31'
scenario="HIST" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,rcode ,tsub_root ,meanVar, str(nclust), outname, str(Ngrid)]
subprocess.check_output(cmd)

# compile timeseries meand RCP26 NEAR
start='2030-01-01' 
end = '2050-12-31'
scenario="RCP26" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,rcode ,tsub_root ,meanVar, str(nclust), outname, str(Ngrid)]
subprocess.check_output(cmd)


# compile timeseries meand RCP26 FAR
start='2080-01-01' 
end = '2099-12-31'
scenario="RCP26" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,rcode ,tsub_root ,meanVar, str(nclust), outname, str(Ngrid)]
subprocess.check_output(cmd)


# compile timeseries meand RCP85 NEAR


start='2030-01-01' 
end = '2050-12-31'
scenario="RCP85" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,rcode ,tsub_root ,meanVar, str(nclust), outname, str(Ngrid)]
subprocess.check_output(cmd)


# compile timeseries meand RCP85 FAR
start='2080-01-01' 
end = '2099-12-31'
scenario="RCP85" 
meanVar= root+"/mean_ts_"+str(col)+"_"+scenario+start+end+".csv"
outname=scenario+"_" +str(col)+"_"+start+"_"+end

timeseries_means(root, nsims, col, start, end, scenario)
cmd = ["Rscript" ,rcode ,tsub_root ,meanVar, str(nclust), outname, str(Ngrid)]
subprocess.check_output(cmd)
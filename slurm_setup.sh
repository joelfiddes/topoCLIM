#!/bin/bash
# JobArray.sh
# run : sbatch slurm_setup.sh /home/caduff/sim/tclim_ch/ /home/caduff/sim/ch_tmapp_50/ /home/caduff/sim/raw_cordex  10

# wd=sys.argv[1] #'/home/joel/sim/qmap/topoclim_ch/'
# raw_dir = sys.argv[2] # /home/caduff/sim/tclim/raw_cordex
# tscale_sim_dir =sys.argv[3] 
# num_cores=sys.argv[4] #10


#SBATCH -J etup # A single job name for the array
#SBATCH -p node # Partition (required)
#SBATCH -A node # Account (required)
#SBATCH -q normal # QOS (required)
#SBATCH -n 10 # one cores
#SBATCH -t 01:00:00 # Running time of 2 days
#SBATCH --mem 4000 # Memory request of 4 GB
#SBATCH -o LOG_setup.out # Standard output - write the console output to the output folder %A= Job ID, %a = task or Step ID
#SBATCH -e LOG_setup.err # Standard error -write errors to the errors folder and
#SBATCH --array=1 # create a array from 1to16 and limit the concurrent runing task  to 50
#SBATCH --mail-user=joelfiddes@gmail.com
#SBATCH --mail-type=ALL  # Send me some mails when jobs end or fail.

pwd; hostname; date

# run sequentially
# $1 is wd
python tclim_hpc_setup.py $1 $2 $3 $4


date



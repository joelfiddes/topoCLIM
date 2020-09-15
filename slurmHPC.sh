#!/bin/sh
#SBATCH -J tmapp
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joelfiddes@gmail.com
#SBATCH --ntasks=50       # tasks requested
#SBATCH --mem-per-cpu=6000
#SBATCH -o outfile  # send stdout to outfile
#SBATCH -e errfile  # send stderr to errfile
#SBATCH -t 30:00:00  # time requested in hour:minute:second

python run_qmap_setup_HPC.py  50
#!/bin/sh
#SBATCH -J tclim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=joelfiddes@gmail.com
#SBATCH --ntasks=5	 # tasks requested
#SBATCH --mem-per-cpu=6000
#SBATCH -o outfile  # send stdout to outfile
#SBATCH -e errfile  # send stderr to errfile
#SBATCH -t 6:00:00  # time requested in hour:minute:second

module load python-3.7.6-gcc-9.1.0-2i2j24b

python3 run_qmap_setup_HPC.py  5

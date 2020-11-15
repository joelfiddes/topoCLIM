#!/bin/bash
# JobArray.sh
# run : sbatch slurm_postqmap.sh  /home/caduff/sim/tclim_ch/ /home/caduff/sim/ch_tmapp_50/  1100
#SBATCH -J postqmap # A single job name for the array
#SBATCH -p node # Partition (required)
#SBATCH -A node # Account (required)
#SBATCH -q normal # QOS (required)
#SBATCH -n 1 # one cores
#SBATCH -t 01:00:00 # Running time of 2 days
#SBATCH --mem 4000 # Memory request of 4 GB
#SBATCH -o LOG_postqmap.out # Standard output - write the console output to the output folder %A= Job ID, %a = task or Step ID
#SBATCH -e LOG_postqmap.err # Standard error -write errors to the errors folder and
#SBATCH --array=1-100 # create a array from 1to16 and limit the concurrent runing task  to 50
#SBATCH --mail-user=joelfiddes@gmail.com
#SBATCH --mail-type=ALL  # Send me some mails when jobs end or fail.

pwd; hostname; date

#Set the number of runs that each SLURM task should do
PER_TASK=$(($3/100))

# Calculate the starting and ending values for this task based
# on the SLURM task and the number of runs per task.
START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

# Print the task and run range
echo This is task $SLURM_ARRAY_TASK_ID, which will do runs $START_NUM to $END_NUM
python tclim_hpc_postqmap.py $1 $2 $START_NUM $END_NUM


date



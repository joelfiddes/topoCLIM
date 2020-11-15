# Example:

# bash tclim_hpc_run.sh /home/caduff/sim/tclim_ch/ /home/caduff/sim/ccamm_inter  /home/caduff/sim/ch_tmapp_50/ /home/caduff/sim/raw_cordex 1100

# Args:
#	$1: tclim workdir
# 	$2: tmapp workdir
#   $3: cordex raw dir
# 	$4: number of samples rounded up to nearest 100

NGRIDS=100

# activate env
source tclim3/bin/activate

# clear logs
rm LOG*

# integer arg is number of cores used by joblib and can be fixed at 10
SBATCHID=$(sbatch slurm_setup.sh $1 $2 10)
jid1=${SBATCHID//[!0-9]/}


SBATCHID=$(sbatch  --dependency=afterany:$jid1 --array=1-$NGRIDS  slurm_qmap.sh $1 $2 $3 $4)
jid2=${SBATCHID//[!0-9]/}

SBATCHID=$(sbatch  --dependency=afterany:$jid2 --array=1-$NGRIDS  slurm_postqmap.sh $1 $2 $4)
jid3=${SBATCHID//[!0-9]/}
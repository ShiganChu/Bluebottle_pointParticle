#!/bin/sh
#
#SBATCH --partition=devel
#SBATCH --gres=gpu:1
#SBATCH --job-name=bb-point
#SBATCH --dependency=afterok:23929
sleep 10
./build-restart.sh $SLURM_JOB_NODELIST $SLURM_JOB_NAME $SLURM_JOB_ID
sbatch submit-restart.sh
srun ./bluebottle -r

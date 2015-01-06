#!/bin/sh

# submit.sh

#

#SBATCH --partition=devel

#SBATCH --gres=gpu:1

#SBATCH --job-name=pointPart

#SBATCH --mail-type=FAIL

#SBATCH --mail-user=chushigan@gmail.com

#SBATCH --output=good

#SBATCH --open-mode=append

srun ./bluebottle >good

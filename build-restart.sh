#!/bin/bash

# build-restart.sh

echo "#!/bin/sh" > submit-restart.sh

echo "#" >> submit-restart.sh

if [ $1 == "guinevere.camelot.me.jhu.edu" ];

then

  echo "#SBATCH --partition=devel" >> submit-restart.sh;

elif [ $1 == "arthur.camelot.me.jhu.edu" ];

then

  echo "#SBATCH --partition=tesla" >> submit-restart.sh;

fi

echo "#SBATCH --gres=gpu:1" >> submit-restart.sh

echo "#SBATCH --job-name=$2" >> submit-restart.sh

echo "#SBATCH --dependency=afterok:$3" >> submit-restart.sh

echo "#SBATCH --output=good" >> submit-restart.sh

echo "#SBATCH --open-mode=append" >> submit-restart.sh

echo "sleep 10" >> submit-restart.sh

echo "./build-restart.sh \$SLURM_JOB_NODELIST \$SLURM_JOB_NAME \$SLURM_JOB_ID" >> submit-restart.sh

echo "sbatch submit-restart.sh" >> submit-restart.sh

echo "srun ./bluebottle -r" >> submit-restart.sh

chmod +x submit-restart.sh

#!/bin/bash 

#SBATCH -p gpu-bf
#SBATCH --gres=gpu:1 
#SBATCH -c 1 
#SBATCH --mem=8g 
#SBATCH -t 0:18:00
#SBATCH -J param_diffuse

#source activate /software/conda/envs/SE3nv/

CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" tasks)
# tell bash to run $CMD
echo "${CMD}"
echo "${CMD}" | bash

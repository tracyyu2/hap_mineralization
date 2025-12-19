#!/bin/bash
#SBATCH -p gpu-bf
#SBATCH --mem=8g
#SBATCH -t 00:20:00
#SBATCH --gres=gpu:1
#SBATCH --exclude=gpu[3-24],gpu[71-83],gpu[110-141],gpu[51-60]
#SBATCH -n 2
echo "Hello from job $SLURM_JOB_ID on $(hostname) at $(date)"
echo $CUDA_VISIBLE_DEVICES
#source activate SE3
CMD=$(head -n $SLURM_ARRAY_TASK_ID tasks | tail -1)
exec ${CMD}

#####command:
#####sbatch -a 1-$(cat tasks|wc -l) run.sh

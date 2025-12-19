#!/bin/bash
#SBATCH -p gpu
#SBATCH --gres=gpu:a4000:1
#SBATCH --mem=16g 
#SBATCH -c 4
#SBATCH -t 14:00:00
#SBATCH -o out

CMD=$(head -n $SLURM_ARRAY_TASK_ID task.list | tail -1)
exec ${CMD}

#sbatch -a 1-$(cat task.list |wc -l) run_gpu.sh

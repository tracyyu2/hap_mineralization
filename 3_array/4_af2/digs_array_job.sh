#!/bin/bash
#SBATCH -J alphafold
#SBATCH -p gpu-bf
#SBATCH --mem=16g
#SBATCH -t 00:10:00
#SBATCH -c 2
#SBATCH --gres=gpu
#SBATCH --exclude=g2513
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tracyy2@uw.edu
echo "Hello from job $SLURM_JOB_ID on $(hostname) at $(date)"
echo $CUDA_VISIBLE_DEVICES
#source activate SE3
CMD=$(head -n $SLURM_ARRAY_TASK_ID tasks | tail -1)
exec ${CMD}

#####command:
#####sbatch -a 1-$(cat tasks|wc -l) run.sh

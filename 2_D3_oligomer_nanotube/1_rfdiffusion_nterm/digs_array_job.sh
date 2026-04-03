#!/bin/bash
#SBATCH -J rfd
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH -t 00:20:00
#SBATCH -c 2
#SBATCH --gres=gpu:a4000:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tracyy2@uw.edu
#SBATCH --array=1-9

# Get the task command from the tasks file
task=$(sed -n "${SLURM_ARRAY_TASK_ID}p" task.list)

# Execute the task command
eval $task

#!/bin/bash
#SBATCH -p cpu
#SBATCH -c 4
#SBATCH -t 10:00
#SBATCH -o out

CMD=$(head -n $SLURM_ARRAY_TASK_ID task.list | tail -1)
exec ${CMD}

#sbatch -a 1-$(cat task.list|wc -l) run_cpu.sh

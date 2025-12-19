#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=2g
#SBATCH -t 01:00:00
/software/containers/pyrosetta.sif realign.py --input_directory ../2_select/08_22_2024_d5_PO4 --seed_directory ../0_seeds/08_22_2024_d5_PO4 --output_directory 08_22_2024_d5_PO4 --task_id $SLURM_ARRAY_TASK_ID
#sbatch -a 1-120 runPos.sub

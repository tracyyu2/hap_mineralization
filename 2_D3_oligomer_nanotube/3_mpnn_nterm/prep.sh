#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=2g
#SBATCH -t 4:00:00
#SBATCH -c 2
#SBATCH --output=prep_run.out
/software/containers/pyrosetta.sif parse_pdbs.py --input_directory ../02_select_nterm/12_12_2024_run1_nterm --output_directory 12_12_2024_run1_nterm

#!/bin/bash
#SBATCH -p cpu
#SBATCH --mem=2g
#SBATCH -t 4:00:00
#SBATCH -c 2
#SBATCH --output=prep_run.out
/software/containers/pyrosetta.sif parse_pdbs.py --input_directory ../02_inspect_diffusion/20251001_rempnn_array53 --output_directory 20251001_rempnn

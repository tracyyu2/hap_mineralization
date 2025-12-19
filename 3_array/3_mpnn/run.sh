#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH -t 00:30:00
#SBATCH --gres=gpu:a4000:1
#SBATCH -c 3
#SBATCH --output=mpnn_run.out

source activate mlfold
export PATH=/home/tracyy2/mpnn-master/:$PATH
/software/containers/mlfold.sif /home/tracyy2/mpnn-master/mpnn_run_tied.py \
        --max_length 10000 \
        --checkpoint_path '/projects/ml/struc2seq/data_for_complexes/training_scripts/paper_experiments/model_outputs/p10/checkpoints/epoch51_step255000.pt' \
        --hidden_dim 192 \
        --num_layers 3 \
        --protein_features 'full' \
        --jsonl_path='pdbs.jsonl' \
        --chain_id_jsonl 'pdbs_masked.jsonl'  \
        --fixed_positions_jsonl 'pdbs_fixed.jsonl' \
        --out_folder='output' \
        --num_seq_per_target 24 \
        --sampling_temp="0.2" \
        --batch_size 4 \
        --omit_AAs 'MC' \
        --backbone_noise 0.05 \
        --decoding_order 'random' \
        --bias_AA_jsonl 'H:-0.9,M:-2.0' \
        --num_connections 64 \
        --tied_positions_jsonl 'pdbs_tied.jsonl'

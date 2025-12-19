from glob import glob
import argparse
import json
import os
import numpy as np

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input_directory",
    default='./',
    help="Path to MPNN outputs",
)

parser.add_argument(
    "--output_directory",
    default='./',
    help="Path to write fold outputs",
)

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

with open(os.path.join(args.output_directory, 'args.txt'), 'w') as f:
    json.dump(vars(args), f, indent=2)

nseqs = 24
cnt = 0

fout_tasks = open(os.path.join(args.output_directory, 'tasks'), 'w')

for fa in glob(os.path.join(args.input_directory, 'output/alignments/*.fa')):
    with open(fa, 'r') as fin:
        fin.readline()  # skip header
        fin.readline()  # skip second line

        for i in range(nseqs):
            out_ind = int(np.floor(cnt / 1000))

            fasta_dir = os.path.join(args.output_directory, f'fastas{out_ind}')
            if not os.path.exists(fasta_dir):
                os.makedirs(fasta_dir)

            fin.readline()  # sequence header
            seq_line = fin.readline().strip()

            #pname = os.path.basename(fa)[:-3] #extracts basename
            pname = fa.split('/')[-1][:-3] #extracts file name
            fasta_filename = f"{pname}_{i}.fa"
            fout_path = os.path.join(fasta_dir, fasta_filename)

            sequences = seq_line.split('/')[:4] #first 4 seqs
            if len(sequences) < 3:
                print(f"Warning: Less than 3 sequences found in {fa}, skipping.")
                continue

            if sequences[0] != sequences[1]:
                print(f"Sequences do not match in {fa} at {i}!!!")
                break

            out_str = ' '.join(sequences) + '\n'

            with open(fout_path, 'w') as fout:
                fout.write(f">{pname}_{i}\n")
                fout.write(out_str)

            task_command = (
                f"/software/containers/users/ikalvet/mlfold3/mlfold3_01.sif "
                f"/net/software/lab/alphafold3/run_alphafold_custom.py "
                f"--fasta {fout_path} "
                f"--run_data_pipeline=false "
                f"--output_dir {args.output_directory}\n"
            )
            fout_tasks.write(task_command)

            cnt += 1

fout_tasks.close()


from glob import glob
import argparse
import json
import os
import numpy as np

def extract_chain_bcd_sequence(line):
    """Extract sequences for chains B, C, and D, assuming sequences are '/'-separated."""
    sequences = line.strip().split('/')
    if len(sequences) >= 4:  # Ensure B, C, and D exist
        return f"{sequences[1]}:{sequences[2]}:{sequences[3]}"  # Join B, C, D with ":"
    return None

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input_directory",
    default='./',
    help="Path to MPNN outputs",
)

parser.add_argument(
    "--output_directory",
    default='./',
    help="Path to write colabfold outputs",
)

args = parser.parse_args()

os.makedirs(args.output_directory, exist_ok=True)

with open(os.path.join(args.output_directory, 'args.txt'), 'w') as f:
    json.dump(args.__dict__, f, indent=2)

nseqs = 24

fout_tasks = open(os.path.join(args.output_directory, 'tasks'), 'w')
cnt = 0

for fa in glob(os.path.join(args.input_directory, 'output/alignments/*.fa')):
    with open(fa, 'r') as fin:
        # Read first two lines (FASTA header and first sequence line)
        fasta_header = fin.readline().strip()
        sequence_line = fin.readline().strip()

        for i in range(nseqs):
            out_ind = int(np.floor(cnt / 1000))  # Ensure integer folder names
            fasta_dir = os.path.join(args.output_directory, f'fastas{out_ind}')
            os.makedirs(fasta_dir, exist_ok=True)

            # Read alternating lines (sequence descriptions and actual sequences)
            fin.readline()  # Skip description line (e.g., ">seq_name")
            sequence_line = fin.readline().strip()

            pname = os.path.basename(fa)[:-3]

            fasta_file = os.path.join(fasta_dir, f"{pname}_{i}.fa")
            with open(fasta_file, 'w') as fout:
                fout.write(f">{pname}_{i}\n")

                chain_bcd_seq = extract_chain_bcd_sequence(sequence_line)
                if chain_bcd_seq:
                    fout.write(chain_bcd_seq + '\n')
                else:
                    print(f'No chain BCD sequence found for {pname}')
                    continue

            fout_tasks.write(
                f"/software/containers/colabfold.sif colabfold_batch "
                f"--msa-mode single_sequence --pair-mode unpaired_paired "
                f"--model-type alphafold2_multimer_v3 --num-recycle 5 "
                f"--recycle-early-stop-tolerance 0.05 --rank multimer "
                f"--model-order 3 --num-models 1 {fasta_file} outputs{out_ind}\n"
            )

            cnt += 1

fout_tasks.close()

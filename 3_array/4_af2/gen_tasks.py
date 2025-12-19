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
    help="Path to write colabfold outputs",
)

args = parser.parse_args()

if not os.path.exists(args.output_directory):
    os.makedirs(args.output_directory)

with open(args.output_directory + '/args.txt', 'w') as f:
    json.dump(args.__dict__, f, indent=2)

nseqs = 24

fout_tasks = open(args.output_directory + '/tasks', 'w')
cnt = 0

for fa in glob(args.input_directory + '/output/alignments/*fa'):
    with open(fa, 'r') as fin:
        lline = fin.readline()
        lline = fin.readline()

        for i in range(nseqs):
            out_ind = np.floor(cnt / 1000)
            if not os.path.exists(args.output_directory + '/fastas%d' % out_ind):
                os.makedirs(args.output_directory + '/fastas%d' % out_ind)

            lline = fin.readline()
            lline = fin.readline()
            pname = fa.split('/')[-1][:-3]

            fout_path = os.path.join(args.output_directory, 'fastas%d' % out_ind, f"{pname}_{i}.fa")
            with open(fout_path, 'w') as fout:
                fout.write(f">{pname}_{i}\n")

                # Only take the first four sequences separated by '/'
                sequences = lline.split('/')[:4]
                if len(sequences) < 4:
                    print(f"Warning: Less than 4 sequences found in {fa}, skipping.")
                    continue

                out_str = ":".join(sequences)

                if sequences[0] != sequences[1]:
                    print("Sequences do not match!!!")
                    break

                fout.write(out_str)
            
            fout_tasks.write("/software/containers/colabfold.sif colabfold_batch " +
                             '--msa-mode single_sequence --pair-mode unpaired_paired ' +
                             '--model-type alphafold2_multimer_v3 --num-recycle 5 ' +
                             '--recycle-early-stop-tolerance 0.05 --rank multimer ' +
                             "--model-order 3 --num-models 1 fastas%d/" % out_ind + 
                             f"{pname}_{i}.fa outputs%d\n" % out_ind)

            cnt += 1

fout_tasks.close()


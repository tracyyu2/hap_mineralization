input_pdb="$1"
num_designs="$2"
weight_intra="$3"
weight_inter="$4"
len_N="$5"
len_C="$6"
keep_start="$7"
keep_end="$8"
thread="$9"

input_name=$(echo "$input_pdb" | sed s'/.pdb//'g)

outname=""$input_name"_intra"$weight_intra"_inter"$weight_inter"_N"$len_N"_C"$len_C"_"$thread""

mkdir -p outputs
outpath="outputs/"

/software/containers/SE3nv.sif \
/net/databases/diffusion/github_repo/rf_diffusion/run_inference.py \
--config-name=symmetry \
inference.model_runner='NRBStyleSelfCond' \
inference.symmetry=\"D2\" \
inference.num_designs=${num_designs} \
inference.output_prefix="${outpath}/${outname}" \
inference.softmax_T=0.0001 \
diffuser.T=50 \
diffuser.aa_decode_steps=0 \
diffuser.chi_type='interp' \
potentials.guiding_potentials=[\"type:olig_contacts,weight_intra:${weight_intra},weight_inter:${weight_inter}\"] \
potentials.olig_intra_all=True \
potentials.olig_inter_all=True \
potentials.guide_scale=2 \
potentials.guide_decay=\"cubic\" \
contigmap.contigs=[\'"$len_N"-"$len_N",A"$keep_start"-"$keep_end","$len_C"-"$len_C",0\ "$len_N"-"$len_N",B"$keep_start"-"$keep_end","$len_C"-"$len_C",0\ "$len_N"-"$len_N",C"$keep_start"-"$keep_end","$len_C"-"$len_C",0\ "$len_N"-"$len_N",D"$keep_start"-"$keep_end","$len_C"-"$len_C"\'] \
inference.input_pdb=\'${input_pdb}\' \
2> ${outpath}/run.log

: '
originally kept
A"$keep_start"-"$keep_end",${insert_len},A128-372
68 resis

try
A"$keep_start"-"$keep_end",${insert_len},A132-372
72 resis
'

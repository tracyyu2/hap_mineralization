import os

# -----------------------------
# Helper function: get PDB length
# -----------------------------
def get_pdb_length(pdb_file):
    residues = set()
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    resi = int(line[22:26])
                    residues.add(resi)
                except ValueError:
                    continue
    return max(residues) if residues else 0


# -----------------------------
# Define parameters
# -----------------------------
input_directory = "/home/tracyy2/Projects/silica_silicate/two_component_designs/0_homo_d3"

input_pdbs = [
    os.path.join(input_directory, f)
    for f in os.listdir(input_directory)
    if f.endswith(".pdb")
]

num_designs = 1
num_repeats = 3   # or 10, or however many designs you want per condition
weight_intras = 1
weight_inters = 1

len_Ns = [20, 22, 24, 26, 28, 30]   # can be []
len_Cs = [0]

# -----------------------------
# Generate thread indices
# -----------------------------
if len(len_Ns) == 0:
    total = len(input_pdbs) * len(len_Cs) * num_repeats
else:
    total = len(input_pdbs) * len(len_Ns) * len(len_Cs) * num_repeats

threads = list(range(1, total + 1))


# -----------------------------
# Generate tasks
# -----------------------------
tasks = []
thread_idx = 0

for input_pdb in input_pdbs:

    pdb_len = get_pdb_length(input_pdb)
    keep_start = 1
    keep_end = pdb_len

    if len_Ns:
        for len_N in len_Ns:
            for len_C in len_Cs:
                for repeat in range(num_repeats):  #number of jobs
                    thread = threads[thread_idx]
                    tasks.append((
                        input_pdb,
                        1,  # num_designs = 1
                        weight_intras,
                        weight_inters,
                        len_N,
                        len_C,
                        keep_start,
                        keep_end,
                        thread
                    ))
                    thread_idx += 1
    else:
        for len_C in len_Cs:
            for repeat in range(num_repeats):  # ID of jobs
                thread = threads[thread_idx]
                tasks.append((
                    input_pdb,
                    1,
                    weight_intras,
                    weight_inters,
                    None,
                    len_C,
                    keep_start,
                    keep_end,
                    thread
                ))
                thread_idx += 1


# -----------------------------
# Write task file
# -----------------------------
os.makedirs("tasks", exist_ok=True)

with open("tasks/task_list.txt", "w") as f:
    for task in tasks:
        input_pdb, num_design, weight_intras, weight_inters, len_N, len_C, keep_start, keep_end, thread = task

        len_N = len_N if len_N is not None else 0

        f.write(
            f"bash ./d3_diffuse.sh "
            f"{os.path.basename(input_pdb)} "
            f"{num_design} {weight_intras} {weight_inters} "
            f"{len_N} {len_C} "
            f"{keep_start} {keep_end} "
            f"{thread}\n"
        )

print(f"Generated {len(tasks)} tasks in tasks/task.list")

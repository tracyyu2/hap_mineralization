import os

# Define parameters
input_directory = "/home/tracyy2/fiber_design/Dihedral/00_dihedral_rosetta_nterm/"
input_pdbs = [os.path.join(input_directory, f) for f in os.listdir(input_directory) if f.endswith(".pdb")]
num_designs = 5
weight_intras = 1
weight_inters = 1
len_Ns = [25, 30, 35, 40]  # If len_Ns should be empty, leave it as an empty list.
len_Cs = [0]  # Multiple length values for len_C to test
keep_start = 1
keep_end = 122  # (start, end)

# Check if len_Ns is empty and adjust task generation accordingly
if len(len_Ns) == 0:
    # Only generate tasks based on len_Cs
    threads = list(range(1, len(input_pdbs) * len(len_Cs) + 1))
else:
    # Generate threads if len_Ns is not empty
    threads = list(range(1, len(input_pdbs) * len(len_Ns) * len(len_Cs) + 1))

# Generate combinations
tasks = []
thread_idx = 0
for input_pdb in input_pdbs:
    if len_Ns:
        # If len_Ns is not empty, iterate over both len_N and len_C
        for len_N in len_Ns:
            for len_C in len_Cs:
                thread = threads[thread_idx]
                tasks.append((input_pdb, num_designs, weight_intras, weight_inters, len_N, len_C, keep_start, keep_end, thread))
                thread_idx += 1
    else:
        # If len_Ns is empty, only iterate over len_C
        for len_C in len_Cs:
            thread = threads[thread_idx]
            tasks.append((input_pdb, num_designs, weight_intras, weight_inters, None, len_C, keep_start, keep_end, thread))
            thread_idx += 1

# Output tasks to a file with the correct format
os.makedirs("tasks", exist_ok=True)
with open("tasks/task_list.txt", "w") as f:
    for task in tasks:
        input_pdb, num_design, weight_intras, weight_inters, len_N, len_C, keep_start, keep_end, thread = task
        # Replace None with 0 for len_N
        len_N = len_N if len_N is not None else 0
        # Write the task in the desired format
        f.write(f"bash ./d3_diffuse.sh {os.path.basename(input_pdb)} {num_design} {weight_intras} {weight_inters} {len_N} {len_C} {keep_start} {keep_end} {thread}\n")

print(f"Generated {len(tasks)} tasks in tasks/task_list.txt")




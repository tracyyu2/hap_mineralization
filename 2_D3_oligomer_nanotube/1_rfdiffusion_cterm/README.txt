#This step is to use Rosetta symmetric diffusion to generate backbones filling the gap between two C3 rings
#change the input and output directories in gen_tasks.py to create a list of tasks
#sbatch -a 1-$(cat task.list |wc -l) run_gpu.sh

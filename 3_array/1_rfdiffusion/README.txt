#this step is to run rf symmetric diffusion to generate backbones filling the gaps bewteen d2 interface
#change the input and output directories in gen_tasks.py to create a list of tasks
#sbatch -a 1-$(cat task.list |wc -l) run_gpu.sh

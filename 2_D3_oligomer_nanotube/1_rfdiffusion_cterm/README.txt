#change the input and output directories in gen_tasks.py to create a list of tasks
#sbatch -a 1-$(cat task.list |wc -l) run_gpu.sh

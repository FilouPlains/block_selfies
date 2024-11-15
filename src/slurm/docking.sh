#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=02:00:00
#SBATCH --job-name=GNU_docking
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/GNU_docking_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/GNU_docking_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024M
#SBATCH --array=0-1999

### make sure good module are loaded
module purge
module restore Blocks_modules

### keep track of modules and python packages
module list
pip list

parallel --jobs $SLURM_CPUS_PER_TASK --joblog logs/parallel_jobs/tasks_${SLURM_ARRAY_TASK_ID}.log < results/mol_mutations/celecoxib/celecoxib_3D_structure/GNU_instructions/task_${SLURM_ARRAY_TASK_ID}.txt

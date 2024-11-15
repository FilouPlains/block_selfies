#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:05:00
#SBATCH --job-name=gen_3D_struct
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/gen_3D_struct/output_logs/gen_3D_struct_%A_%a.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/gen_3D_struct/error_logs/gen_3D_struct_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024M
#SBATCH --array=0-1999

### make sure good module are loaded
module purge
module restore Blocks_modules

### keep track of modules and python packages
module list
pip list

### skip line for log cleanliness
echo

### launch script
python scripts/gen_3D_struct.py --input_csvpath results/mol_mutations/celecoxib/celecoxib_3D_structure/subset_${SLURM_ARRAY_TASK_ID}.csv 

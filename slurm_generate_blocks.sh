#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:58:00
#SBATCH --job-name=blocks_generation
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/generate_blocks/output_logs/generate_blocks_%A_%a.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/generate_blocks/error_logs/generate_blocks_%A_%a.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024M
#SBATCH --array=0-3999

### make sure good module are loaded
module purge
module restore Blocks_modules

### keep track of modules and python packages
module list
pip list

### skip line for log cleanliness
echo

### launch script
python scripts/smiles2blocks.py data/SMILES_data/subsets_dir/subset_${SLURM_ARRAY_TASK_ID}.csv 4 1000 25

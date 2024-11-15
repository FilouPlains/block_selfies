#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:40:00
#SBATCH --job-name=blocks_generation
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/generate_blocks_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/generate_blocks_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024M


### make sure good module are loaded
module purge
module restore Blocks_modules

### keep track of modules and python packages
module list
pip list

### skip line for log cleanliness
echo

### launch script
python scripts/smiles2blocks.py data/SMILES_data/subsets_dir/subset_500.csv 4 1000 25
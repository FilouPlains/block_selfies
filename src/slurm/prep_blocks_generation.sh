#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:5:00
#SBATCH --job-name=prep_blocks_generation
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/prep_blocks_generation_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/prep_blocks_generation_%A.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=40000M

### make sure good module are loaded 
module purge
module restore Blocks_modules

### keep track of modules and python pacakges 
module list
pip list 

### skip a line for log cleanliness
echo

### launch download script 
python scripts/prep_blocks_generation.py
#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=01:00:00
#SBATCH --job-name=mutate_mol
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/mutate_mol_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/mutate_mol_%A.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=40000M

### make sure good module are loaded 
module purge
module restore Blocks_modules

### keep track of modules and python pacakges 
module list
pip list 

### launch download script 
python scripts/mutate_mol.py 
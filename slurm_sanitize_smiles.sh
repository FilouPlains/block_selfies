#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:20:00
#SBATCH --job-name=sanitize_smiles
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/sanitize_smiles_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/sanitize_smiles_%A.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=10000M

### make sure good module are loaded 
module purge
module restore Blocks_modules

### keep track of modules and python pacakges 
module list
pip list 

### skip a line for log cleanliness
echo

### launch  script 
python scripts/sanitize_smiles.py

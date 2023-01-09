#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:10:00
#SBATCH --job-name=download_MOSES
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/download_MOSES_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/download_MOSES_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000M

### make sure good module are loaded 
module purge
module restore whole_canonical_modules

### keep track of modules and python pacakges 
module list
pip list 

### skip a line for log cleanliness 
echo

### launch download script 
python scripts/download_MOSES.py 
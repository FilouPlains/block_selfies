#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=01:00:00
#SBATCH --job-name=concatenate2lib
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/concatenate2lib_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/concatenate2lib_%A.err
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

### launch  script 
python scripts/concatenate2lib.py

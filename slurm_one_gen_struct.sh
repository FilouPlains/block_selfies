#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:30:00
#SBATCH --job-name=blocks_generation
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/gen_3D_struct/output_logs/gen_3D_struct_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/gen_3D_struct/error_logs/gen_3D_struct_%A.err
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
python scripts/gen_3D_struct.py results/mol_mutations/celecoxib/celecoxib_3D_structure/subset_0.csv 

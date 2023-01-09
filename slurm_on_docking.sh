#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:20:00
#SBATCH --job-name=GNU_docking
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/GNU_docking_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/GNU_docking_%A.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=1024M

### make sure good module are loaded
module purge
module restore Blocks_modules

### keep track of modules and python packages
module list
pip list

python scripts/docking_pipeline.py --ligand_path results/mol_mutations/celecoxib/celecoxib_3D_structure/structures_files/mol_42.pdbqt
# parallel --jobs 4 --joblog logs/parallel_jobs/tasks_0.log < results/mol_mutations/celecoxib/celecoxib_3D_structure/GNU_instructions/task_0.txt

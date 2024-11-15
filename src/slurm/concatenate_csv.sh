#!/bin/sh
#SBATCH --account=def-jeromew
#SBATCH --time=00:40:00
#SBATCH --job-name=rename_file
#SBATCH --output=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/rename_file_%A.out
#SBATCH --error=/home/retienne/projects/def-jeromew/retienne/blocks_SELFIES/logs/rename_file_%A.err
#SBATCH --cpus-per-task=20
#SBATCH --mem=10000M

python scripts/concatenate_csv.py
#!/bin/bash -l
#SBATCH --job-name=get_ML_data
#SBATCH --output=data_init.txt
#SBATCH --error=data_init_err.txt
#SBATCH --time=01:00:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --account="project_xxxxxxx"

module load geoconda

python -u get_data.py

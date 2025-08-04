#!/bin/bash -l
#SBATCH --job-name=run_ML_basis
#SBATCH --output=ML_output_%j.txt
#SBATCH --error=ML_error_%j.txt
#SBATCH --time=06:10:00
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --account="project_xxxxxx"

set -euo pipefail

module load geoconda

# Run your Python script with Dask

python ML_basis.py

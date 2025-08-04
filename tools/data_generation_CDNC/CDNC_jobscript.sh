#!/bin/bash
#SBATCH --job-name=CDNC_multi_job_submit
#SBATCH --output=CDNC_multi_job_submit_%j.out
#SBATCH --error=CDNC_multi_job_submit_%j.err
#SBATCH --account="project_xxxxxx"              # Name of the project
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:01:00

# array of months (as integers in range from 1 to 12)
MONTHS=(1 2 3 4 5 6 7 8 9 10 11 12)

# submit each month as individual job
for M in "${MONTHS[@]}"; do
    sbatch CDNC_runner.sh $M
done



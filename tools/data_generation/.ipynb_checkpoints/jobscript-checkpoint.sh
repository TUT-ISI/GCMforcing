#!/bin/bash
#SBATCH --job-name=multi_job_submit
#SBATCH --account="project_2010692"              # Name of the project
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:00:30

# array of months (as integers in range from 1 to 12)
# old (1 2 3 4 5 6 7 8 9 10 11 12)
MONTHS=(1)

# submit each month as individual job
for M in "${MONTHS[@]}"; do
    sbatch the_runner.sh $M
done


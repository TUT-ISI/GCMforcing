#!/bin/bash -l
#SBATCH --job-name="3GlobalTraj-CE_CDNC"              # Name of the job
#SBATCH --output=slurm_GlobalTraj-CE_CDNC_%j.txt      # Outputfile
#SBATCH --error=slurm_GlobalTraj-CE_CDNC_%j.txt       # Error file
#SBATCH --account="project_xxxxxx"              # Name of the project
#SBATCH --partition=small                        # Partition the job is submitted for: 2 or more tasks require large partition
#SBATCH --time=01:00:00                          # Maximum runtime (change according to number of parameter combinations assigned in the_runner.py-file)
#SBATCH --ntasks=1                               # Number of parallel jobs
#SBATCH --cpus-per-task=6                       # Number of CPU cores per job
#SBATCH --mem-per-cpu=6G                         # Memory reserved for one core

# # Load necessary modules
# module load python-data
# module load cdo
# module load geoconda

# set current directory (the one containing the_runner.py and other files needed for running it)
cd ./tools/data_generation_CDNC

# number of month for this iteration (from jobscript)
MONTH=$1

# Run the calculations with different sets of parameters
srun -n 1 CDNC_driver.sh $MONTH

# save memory usage at the end of the job to a file
seff $SLURM_JOBID > ${MONTH}memory_usage.txt

#!/bin/bash -l
#SBATCH --job-name="3GlobalTraj-CE"              # Name of the job
#SBATCH --output=slurm_GlobalTraj-CE_%j.txt      # Outputfile
#SBATCH --error=slurm_GlobalTraj-CE_%j.txt       # Error file
#SBATCH --account="project_2010692"              # Name of the project
#SBATCH --partition=small                        # Partition the job is submitted for: 2 or more tasks require large partition
#SBATCH --time=00:15:00                          # Maximum runtime (change according to number of parameter combinations assigned in the_runner.py-file)
#SBATCH --ntasks=1                               # Number of parallel jobs
#SBATCH --cpus-per-task=6                        # Number of CPU cores per job
#SBATCH --mem-per-cpu=6G                         # Memory reserved for one core
#SBATCH --mail-type=ALL                          # Send email at job completion
#SBATCH --mail-user=atte.laakso@uef.fi              # Email address for notifications

# Set current directory (the one containing the_runner.py and other files needed for running it)
cd /scratch/project_2010692/tools/data_generation

# Number of month for this iteration (from jobscript)
MONTH=$1

# Run the calculations with different sets of parameters
srun -n 1 driver.sh $MONTH

# Save memory usage at the end of the job to a file
seff $SLURM_JOBID > ${MONTH}memory_usage.txt

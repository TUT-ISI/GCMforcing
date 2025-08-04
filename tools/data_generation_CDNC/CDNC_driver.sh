#!/bin/bash
# Load necessary modules
module load cdo
module load geoconda

# Read the arguments passed to driver.sh
month=$1

# Pass the arguments to the Python script
python3 CDNC_runner.py $month

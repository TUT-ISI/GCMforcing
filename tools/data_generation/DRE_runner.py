#!/usr/bin/python3
"""
by  Atte Laakso / Aalto University

For running DRE_offline_driver based on parameters defined in this file

NOTE:
this program is designed for running calculations on compute nodes on CSC's Puhti.
"""

# +
import sys, os
import signal
import numpy as np
from functools import partial
from itertools import islice
from itertools import product
from multiprocessing import Pool, freeze_support, cpu_count
from write_netcdf import *
from DRE_offline_driver import running
from DRE_combinations import make_combos
import logging
# initialize logging for later debugging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# Models for RH: "CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","INCA","GISS-ModelE2p1p1-OMA","MIROC-SPRINTARS","OsloCTM3v1.01"
# Models with mmrC: "CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA","GFDL-AM4"
    
# Setting up constants for the offline_driver
light = "550"     # Wavelenght of light in nanometers (550 by standard)

"""
Define the models and parameters to be used
"""
# Base model that other models are compared to
# (check that these are stored in the right folder spesified in the offline_driver)
model1 = "vbs_sensitivity_base_2010.01" 
    
# Comparable models (name these to match the formatting given in the documentation)
# (ensure that every model in this list have all the files corresponding to each parameters stated below (for refractive indices that the values are in the table))
# all current models: ["CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA"]
models = ["CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA","GFDL-AM4"]

# Timeline for sbatch script to pick month
# in format "01","02","03","04","05","06","07","08","09","10","11","12"
time = ["01","02","03","04","05","06","07","08","09","10","11","12"]

# For allowing changes in the rh field to the models,type 'rh', otherwise type 'norh'
rhsel = 'rh'

# For allowing changes in the mmr field to the models, type 'mmr', otherwise type 'nommr'
mmrsel = 'mmr'

# For allowing changes in the refractive indices to the ones of the model, type 'refrac', otherwise type 'nor'
refsel = 'refrac'

# For allowing changes in the total cloud fraction to the ones of the model, type 'clt', otherwise type 'noclt'
cltsel = 'clt'

# For avoiding unneccessary iterations, select True in order to block the code for recalculating the conditions that are already calculated
avoid_print_on = False

# Number of CPUs given to this task (set the same value as in the_runner.sh cpu-per-task)
num_cpus = 6

# Path to Output folder
outPath = '../tools/output_DRE'
# Path to SALSA base data folder
SALSAPath = '../SALSA_2010'
# Path to temporary files
temp = '../tools/temporary'
# Path to AeroCom models
AeroCom = '../tools/AeroCom_models'

# Run the calculations based on values given above as well as in the SLURM script the_runner.sh
def main(month):
    # Update logging
    logging.info(f'Starting computations for month {month}')
    # List all parameters
    parameters = make_combos(rhsel,mmrsel,refsel,cltsel)
    
    # Defining an argument list for parallel processing
    arg_list = []
    for paras in parameters:
        for m2 in models:
                arg_list.append((model1,m2,time[month-1],light,paras[0],paras[1],paras[2],paras[3],avoid_print_on,[outPath,SALSAPath,temp,AeroCom]))
    # Update logging
    logging.info(f'Argument list length: {len(arg_list)}')
    
    # Creating pool of worker processes in order to process parallel given lists of arguments
    with Pool(processes = num_cpus) as pool:
        # Devide work to processable chunks
        stepper=1
        for chunk in chunk_iterations(arg_list, num_cpus):
            # Process one chunk
            logging.info(f'Starting process number {stepper}')
            stepper+=1
            # Using starmapping to call running-method of AOD_offline_driver
            pool.starmap(running, chunk)
    logging.info('Iterations finished')
    
    # Close the pool to new tasks
    pool.close()

    # Wait for all worker processes to complete
    pool.join()



# Chunk the iterations to pieses that are prosessable with amount of parallel processes available to avoid exeeding memory limit
def chunk_iterations(args_list, size):
    # args_list into iterator 
    it = iter(args_list)
    # Creates chunks of the size determined and returns them one-by-one
    while True:
        # Form a chunk
        chunk = list(islice(it, size))
        if not chunk:
            break
        # Return a chunk when needed
        yield chunk


if __name__ == "__main__":
    
    # Month to process
    month=int(sys.argv[1])
    
    # Run main()
    main(month)
# -




#!/usr/bin/python3
import sys, os
import numpy as np
from itertools import islice
from multiprocessing import Pool, freeze_support, cpu_count
from write_netcdf import *
from CDNC_offline_driver import running
from combinations import make_CDNC_combos
import logging
"""
by Atte Laakso / Aalto University

For running CDNC offline driver on compute nodes

NOTE: Updraft has never been read from AeroCom models due to lack of data so far but the functionality is there
"""

# initialize logging for later debugging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

# currently compatable models: 
# "CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA"
# AND "GFDL-AM4" is missing clw file, otherwise would be good to go

# setting up constants for the driver

light = "550"     # wavelenght of light in nanometers (550 by standard)
time = ["01","02","03","04","05","06","07","08","09","10","11","12"] # time values for the system in correct format

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
"""
Define the models and general settings
"""
# base model that other models are compared to
# (check that these are stored in the right folder which is spesified in the CDNC_offline_driver)
model1 = "vbs_sensitivity_base_2010.01" 
    
# comparable models (name these to match the formatting given in the documentation (Naming from Aerocom AP3-CTRL protocol))
# (ensure that every model in this list have all the files corresponding to each of the parameters stated below)
# ["CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA"]
models = ["CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA"]

# for avoiding unneccessary iterations, select True in order to block the code for recalculating the conditions that are already calculated
avoid_print_on = False

# for calculating energy balance with aerosol cloud radiative effect, select true
calculate_CRE = True

# number of CPUs given to this task (when running on compute nodes, set the same value as 
# in CDNC_runner.sh cpu-per-task)
num_cpus = 6

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
"""
Select months and parameters

Note: These values are set different based on whether running on jupyter or on compute nodes
"""
# Define if changes in these parameters are alowed or not (the program creates combinations based on permissions to change certain field)
# To get all the possible parameter combinations, every selection should be positive (eg, 'mmr','nmr')

# for allowing changes in/changing the mmr field to the models, type 'mmr', otherwise type 'nommr'
mmrsel = 'mmr'

# for allowing changes in/changing the nmr field to the models (that is calculated from the mmr data), type 'nmr', otherwise type 'nonmr'
nmrsel = 'nmr'

# for allowing changes in/changing updraft field to the models, type 'updraft', otherwise type 'noupdraft'
updraftsel = 'noupdraft'

# for allowing changes in/changing cloud area fraction field to the models, type 'clt', otherwise type 'noclt'
cltsel = 'clt'

# for allowing changes in/changing cloud liquid water content field to the models, type 'clw', otherwise type 'noclw'
clwsel = 'clw'

# Path to Output folder
outPath = '../tools/output_CDNC'
# Path to SALSA base data folder
SALSAPath = '../SALSA_2010'
# Path to temporary files
temp = '../tools/temporary'
# Path to AeroCom models
AeroCom = '../tools/AeroCom_models'

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

# run the calculations based on values given above or in the SLURM script CDNC_runner.sh
def main(month_sel):
    # update logging
    logging.info(f'Starting computations for month {month_sel} on {num_cpus} cores')
    # list all parameters
    parameters = make_CDNC_combos(mmrsel,nmrsel,updraftsel,cltsel,clwsel)
    
    # defining an argument list for parallel processing
    arg_list = []
    for paras in parameters:
        for m2 in models:
                arg_list.append((model1,m2,time[month_sel-1],light,paras[0],paras[1],paras[2],paras[3],paras[4],avoid_print_on,calculate_CRE,[outPath,SALSAPath,temp,AeroCom]))
    # update logging
    logging.info(f'Argument list length: {len(arg_list)}')
    
    # creating pool of worker processes in order to process parallel given lists of arguments
    with Pool(processes = num_cpus) as pool:
        # devide work to processable chunks
        stepper=1
        for chunk in chunk_iterations(arg_list, num_cpus):
            # process one chunk
            logging.info(f'Starting process number {stepper}')
            stepper+=1
            # using starmapping to call running-method of CDNC_offline_driver
            pool.starmap(running, chunk)
    logging.info('Iterations finished')
    
    # close the pool to new tasks
    pool.close()

    # wait for all worker processes to complete
    pool.join()


# chunk the iterations to pieses that are prosessable with amount of parallel processes available to avoid exeeding memory limit
def chunk_iterations(args_list, size):
    # args_list into iterator 
    it = iter(args_list)
    # creates chunks of the size determined and returns them one-by-one
    while True:
        # form a chunk
        chunk = list(islice(it, size))
        if not chunk:
            break
        # return a chunk when needed
        yield chunk


if __name__ == "__main__":
    # month to process
    month=int(sys.argv[1])
    
    # run main()
    main(month)
# -




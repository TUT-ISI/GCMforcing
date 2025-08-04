#!/usr/bin/python3
#    Atte Laakso
"""
For running different offline drivers

NOTE:
This is the simpler version of the_runner.py made for computing small sets of parameter combinations. 
For larger sets, use the_runner.py via jobscript.sh (on compute nodes)

Please check all the variable definitions above main function according to the suggestions given in the documentation of this code
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
from AOD_offline_driver import running
from offline_driver_modelCDNC import CDNCrunner
from combinations import make_combos
    
# models for RH: "CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","INCA","GISS-ModelE2p1p1-OMA","MIROC-SPRINTARS","OsloCTM3v1.01"
# models with mmr: "CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA","GFDL-AM4"
    

#Setting up constants and testing parameters for the offline_driver
light = "550"     # wavelenght of light in nanometers (550 by standard)
testmodels = ["ECHAM6.3-HAM2.3"]     # name one model to be used in testing
testtime = ["01","02"]     # couple of months to be used in testing


"""
Define the models and parameters to be used
"""

# base model that other models are compared to
# (check that these are stored in the right folder spesified in the offline_driver)
model1 = "vbs_sensitivity_base_2010.01" 
    
# comparable models (name these to match the formatting given in the documentation)
# (ensure that every model in this list have all the files corresponding to each parameters stated below (for refractive indices that the values are in the table))
models = ["CAM5-ATRAS"]

# Select the months for the comparison
# in format "01","02","03","04","05","06","07","08","09","10","11","12"
time = ["01"]

# for changing the rh field to the models,type 'rh', otherwise type 'norh'
rhsel = 'rh'

# for changing the mmr field to the models, type 'mmr', otherwise type 'nommr'
mmrsel = 'nommr'

# for changing the nmr field to the models (that is calculated from the mmr data), type 'nmr', otherwise type 'nonmr'
nmrsel = 'nonmr'

# for changing the refractive indices to the ones of the model, type 'refrac', otherwise type 'nor'
refsel = 'refrac'

#for combining aods for the final output, write True, otherwise False
com_aod = True

# for avoiding unneccessary iterations, select True in order to block the code for recalculating the conditions that are already calculated
avoid_print_on = False

# Get the number of CPUs allocated by SLURM (if wanted to be used)
num_cpus = 2

# define main
def main():
    #defining an argument list
    arg_list = []
    # create list of arguments (all the different parameter combinations)
    for m2 in models:
        for t in time:
            arg_list.append((model1,m2,t,light,rhsel,mmrsel,nmrsel,refsel,com_aod,avoid_print_on))
    
    # creating pool of worker processes in order to process parallel given lists of arguments
    with Pool(processes = num_cpus) as pool:
        #devide work to processable chunks
        stepper=1
        for chunk in chunk_iterations(arg_list, num_cpus):
            print(f'Process number {stepper} has started.')
            stepper+=1
            # process each chunk in the pool
            pool.starmap(running, chunk)
    print("iterations finished")
    
    # close the pool to new tasks
    pool.close()

    # wait for all worker processes to complete
    pool.join()

# for m2 in testmodels:
#     for t in testtime:
#         CDNCrunner(model1,m2,t)
#         print('CDNC iteration '+m2+' '+t+' is done!')



# chunk the iterations to pieses that are prosessable with amount of parallel processes available to avoid exeeding memory limit
def chunk_iterations(args_list, size):
    # args_list into iterator 
    it = iter(args_list)
    # creates chunks of the size determined and returns them one-by-one
    while True:
        chunk = list(islice(it, size))
        if not chunk:
            break
        yield chunk


if __name__ == "__main__":
    main()
# -


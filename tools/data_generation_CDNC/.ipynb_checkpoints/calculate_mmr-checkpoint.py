# Written for Python v3.0
# by Atte Laakso
"""

Calculates the mass mixing ratios' ratio between "comparison model" and SALSA
 
"""

from netCDF4 import Dataset
import sys, os
import numpy as np
import netCDF4 as nc
from write_netcdf import *
from cdo import Cdo as CDO

    #creating method that can be called elsewhere
    #fname_mmrham should involve the model that we want to compare to the basemodel in the file fname_mmrcomp
    #by standard the basemodel should be SALSA
def c_mmr(fname_mmrham,fname_mmrcomp):
    
    # get PM1 for model levels
    # print('Reading in PM1 from file')
    hamdata=nc.Dataset(fname_mmrham, 'r', format='NETCDF4_CLASSIC')    
    pm1 = hamdata.variables['mmrpm1']
    mmrpm1 = np.array(pm1)

    #Same for comparison model
    compdata=nc.Dataset(fname_mmrcomp, 'r', format='NETCDF4_CLASSIC')    
    pm1_comp = compdata.variables['mmrpm1']
    mmrpm1_comp = np.array(pm1_comp)

    #calculating ratio between PM1 in given model and comparison model
    ratio=mmrpm1/mmrpm1_comp

    # close the files
    try:
        hamdata.close()
        compdata.close()
    except IOError as e:
        pass
    
    #returns the ratio
    return ratio
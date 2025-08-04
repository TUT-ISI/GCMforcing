# Written for Python 3
# by Atte Laakso
"""

Calculates the number mixing ratio using mass mixing ratio
 
"""

from netCDF4 import Dataset
import sys, os
import numpy as np
import netCDF4 as nc
from salsa_parameters import *
import math

# Creating method that can be called elsewhere
# fname_mmrham should involve the model that we want to compare to the basemodel in the file fname_mmrcomp
def c_nratio(fname_mmrham,salsa_ham):
    
    # Get PM1 for model levels
    print('Reading in PM1')
    hamdata=nc.Dataset(fname_mmrham, 'r', format='NETCDF4_CLASSIC')    
    pm1 = hamdata.variables['mmrpm1']
    mmrpm1 = np.array(pm1)

    # Close the file and clear unused vriables
    try:
        hamdata.close()
    except IOError as e:
        pass
    del pm1
    
    # Empty dictionary for nmr values
    number_mixing_ratio=dict()


    print('Calculating nmrs')
    # Iterate over all the bins
    for b in bins:
        # Read data from the ham file
        bindata=nc.Dataset(salsa_ham, 'r', format='NETCDF4_CLASSIC')

        # Get the radius of the particles for the bin
        radius=bindata.variables['rdry_'+b]
        radius=np.array(radius)
        
        # Calculate volume of one particle
        aerosol_volumes=4/3*math.pi*(radius**3)
        
        # Read the density of the bin
        bindensity_netc=bindata.variables['densaer_'+b]
        bindensity=np.array(bindensity_netc)

        # Free memory by closing files and deleting variables
        try:
            bindata.close()
        except IOError as e:
            pass
        del radius
        
        # Calculate mass of one aerosol particle
        aerosol_mass=aerosol_volumes*bindensity

        # Calculate the number mixing ratio
        nmr=mmrpm1/aerosol_mass

        # Save nmr to the dicitonary under NUM-variables (for each bin)
        number_mixing_ratio['NUM_'+b]=nmr
        
    #returns the number mixing ratio
    return number_mixing_ratio



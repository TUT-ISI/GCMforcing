"""
Written for Python 3
by Atte Laakso / UEF

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
def c_nratio(fname_mmrham,salsa_ham,densityfile):
    
    # Get PM1 for model levels
    with nc.Dataset(fname_mmrham, 'r', format='NETCDF4_CLASSIC') as hamdata:
        pm1 = hamdata.variables['mmrpm1']
        mmrpm1 = np.array(pm1)

    # Get dry bin densities from file calculated by calculate_dry_density.py (time for this selected in offline driver)
    dry_data = nc.Dataset(densityfile, 'r', format='NETCDF4_CLASSIC')

    # Close the file and clear unused variables
    try:
        hamdata.close()
    except IOError as e:
        pass
    del pm1
    
    # Empty dictionary for nmr values
    number_mixing_ratio=dict()

    # Read data from the ham file
    bindata=nc.Dataset(salsa_ham, 'r', format='NETCDF4_CLASSIC')

    print('Calculating nmrs')
    # Iterate over all the bins
    for b in bins:
        if b in dry_data.variables:
            # Get the radius of the particles for the bin
            radius=bindata.variables['rdry_'+b]
            radius=np.array(radius)
            
            # Calculate volume of one particle
            aerosol_volumes=4/3*math.pi*(radius**3)
            
            # Read the dry aerosol density [kg of aerosol in m³ of aerosol] of the bin
            bindensity_netc = dry_data.variables[b]
            bindensity = np.array(bindensity_netc) # Full read
            print("bindensity shape after read:", bindensity.shape)

            # Calculate mass of one aerosol particle
            aerosol_mass=aerosol_volumes*bindensity

            # Calculate the number mixing ratio
            # Avoid division by zero
            aerosol_mass[aerosol_mass < 1e-30] = np.nan
            nmr = mmrpm1 / aerosol_mass

            # NaN values should be zero
            nmr = np.nan_to_num(nmr)
            nmr = np.clip(nmr, 0, None) # Limit to zero (can't be negative)

            # Save nmr to the dicitonary under NUM-variables (for each bin)
            number_mixing_ratio['NUM_'+b]=nmr

    #returns the number mixing ratio
    return number_mixing_ratio



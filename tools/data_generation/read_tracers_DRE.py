from netCDF4 import Dataset
from parameters import *
from salsa_parameters import *
import numpy as np
import sys, os

def read_aero_bins_dp(tracefilename, vphysfilename, mratio, mmrsel):
    """
    Written for Python v2.7.6
    by  Harri Kokkola / FMI 
        Atte Laakso / UEFs

    Reads in aerosol number mixing ration
    and mass mixing ratio of chemical compounds
    and scales them based on ratio between mass mixing ratios of SALSA and comparison model
    """
    
    # File for reading mixing ratios of individual compounds
    binfile = Dataset(tracefilename, 'r')
    cofile = Dataset(vphysfilename, 'r')

    # Read air density from vphysc file
    rhoam1 = cofile.variables['rhoam1'][:,:,:,:]
    airdensity = rhoam1

    # Volume concentrations
    pvols=dict()
    # Number concentrations
    pnaero=dict()
    pvols_temp=dict()

    for var in variables:
        if var in binfile.variables:
            # Name of tracer species / number concentration
            spec=var[:-4]
            if spec=='NUM':
                if mmrsel == 'mmr':
                    # calculate the number concentration [1/m3] from number mixing ratio of the model
                    pnaero[var]=binfile.variables[var]*airdensity*mratio
                    # limit value for empty cells            
                    pnaero[var]=np.maximum(0.01,pnaero[var])
                else:
                    # calculate the number concentration [1/m3] from number mixing ratio of the SALSA base model (case where only mmr is taken into account) 
                    pnaero[var]=binfile.variables[var]*airdensity
                    # limit value for empty cells            
                    pnaero[var]=np.maximum(0.01,pnaero[var])
            else:
                if mmrsel == 'mmr':
                # calculate volume concentration [m3/m3] from mixing ratios
                    pvols[var]=np.maximum(0.0,binfile.variables[var]*airdensity*mratio/density[spec])
                else:
                    pvols[var]=np.maximum(0.0,binfile.variables[var]*airdensity/density[spec])
        else:
            pvols[var]=airdensity*0.0

    return pnaero,pvols

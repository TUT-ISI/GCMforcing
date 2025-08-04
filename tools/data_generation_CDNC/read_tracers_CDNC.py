# Written for Python v2.7.6
#    Harri Kokkola / FMI 
#
# Reads in aerosol number mixing ration
# and mass mixing ratio of chemical compounds

from netCDF4 import Dataset
from parameters import *
from salsa_parameters import *
import numpy as np
import sys, os

def read_aero_binsdp(tracefilename, vphysfilename, mratio, mmrsel, nmrsel, fname_nmr):
    """
    Reads netCDF file for the (input) bins with the particle concentration data.
    """
    
    # File for reading mixing ratios of individual compounds
    binfile = Dataset(tracefilename, 'r')
    cofile = Dataset(vphysfilename, 'r')

    # time = binfile.variables['time_bnds'][:]

    # Read air density from vphysc file
    rhoam1 = cofile.variables['rhoam1'][:,:,:,:]
    airdensity = rhoam1

    # calculate nmr
    if nmrsel == 'nmr':
        if os.path.exists(fname_nmr):
            with Dataset(fname_nmr, 'r') as nf:
                number_mixing_ratio = dict()
                for b in bins:
                    var_name = f'NUM_{b}'
                    if var_name in nf.variables:
                        number_mixing_ratio[var_name] = np.array(nf.variables[var_name])
                    else:
                        print(f"Warning: {var_name} not found in {fname_nmr}")
        else:
            print("NMR files must be calculated first with calculate_nmr_files.py script!")
            raise FileNotFoundError(f"NMR file not found: {fname_nmr}")


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
                if nmrsel == 'nmr':
                    # calculate the number concentration [1/m3] from number mixing ratio of the model
                    pnaero[var]=number_mixing_ratio[var]*airdensity
                    # limit value for empty cells            
                    pnaero[var]=np.maximum(0.01,pnaero[var])
                else:
                    # calculate the number concentration [1/m3] from number mixing ratio of the SALSA base model (case where only mmr is taken into account) 
                    pnaero[var]=binfile.variables[var]*airdensity*mratio
                    # limit value for empty cells            
                    pnaero[var]=np.maximum(0.01,pnaero[var])
            else:
                if mmrsel == 'mmr':
                # calculate volume concentration [m3/m3] from mixing ratios
                    pvols[var]=np.maximum(0.0,binfile.variables[var]*airdensity*mratio/density[spec])
                else:
                    pvols[var]=np.maximum(0.0,binfile.variables[var]*airdensity/density[spec])
                #if spec=='OC':
                #   pvols_temp[var]=pvols[var]*5.
                #   pvols[var]=pvols_temp[var]
                #   print 'OC found'
        else:

            pvols[var]=airdensity*0.0

    return pnaero,pvols

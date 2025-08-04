# Written for Python v2.7.6
#    Harri Kokkola / FMI 
#
# Reads in aerosol number mixing ration
# and mass mixing ratio of chemical compounds

from netCDF4 import Dataset
from parameters import *
from salsa_parameters import *
import numpy as np
from calculate_nratio import c_nratio

def read_aero_binsdp(tracefilename, vphysfilename, fname_mmrham, fname_ham, lon, lat, mratio, mmrsel, nmrsel):
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
       number_mixing_ratio=c_nratio(fname_mmrham,fname_ham)
    
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

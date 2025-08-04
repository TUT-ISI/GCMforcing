# Written for Python v2.7.6
#    Harri Kokkola / FMI 
#
# Reads in aerosol number mixing ration
# and mass mixing ratio of chemical compounds

from netCDF4 import Dataset
from parameters import *
from salsa_parameters import *
import numpy as np

def read_aero_bins(tracefilename, vphysfilename):
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

   # volume concentrations
   pvols=dict()
   # number concentrations
   pnaero=dict()
   pvols_temp=dict()
   for var in variables:

      if var in binfile.variables:
         # Name of tracer species / number concentration
         spec=var[:-4]
         if spec=='NUM':
            # calculate the number concentration [1/m3] from number mixing ratios 
            pnaero[var]=binfile.variables[var]*airdensity
            # limit value for empty cells            
            pnaero[var]=np.maximum(0.01,pnaero[var])

         else:
            # calculate volume concentration [m3/m3] from mixing ratios
            pvols[var]=np.maximum(0.0,binfile.variables[var]*airdensity/density[spec])
            #if spec=='OC':
            #   pvols_temp[var]=pvols[var]*5.
            #   pvols[var]=pvols_temp[var]
            #   print 'OC found'
      else:

         pvols[var]=airdensity*0.0

   return pnaero,pvols

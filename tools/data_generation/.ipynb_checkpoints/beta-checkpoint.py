from netCDF4 import Dataset
from salsa_parameters import *
from asy_properties import asymmetry
import numpy as np
import sys, os
from cdo import Cdo as CDO
from write_netcdf import *

###############################################################################################
# Code by Atte Laakso 
# based on original program in aod.py by Harri Kokkola

# [1] Taufiq Hassan, H. Moosmüller, Chul E. Chung,
# Coefficients of an analytical aerosol forcing equation determined with a Monte-Carlo radiation model,
# Journal of Quantitative Spectroscopy and Radiative Transfer,
# Volume 164, 2015, Pages 129-136, ISSN 0022-4073,
# https://doi.org/10.1016/j.jqsrt.2015.05.015.
# (https://www.sciencedirect.com/science/article/pii/S0022407315002009)
###############################################################################################


def calculate_beta(nr,ni,dwet,pnaero,pgrheight,lutfile,wavelength):
    """
    Calculates weighted average backscattering fraction beta for given volume concentration fields
    """

    # Open netcdf look-up-table of optical properties
    lookup=Dataset(lutfile, 'r')

    # Mie size parameter
    size_param=dict()

    # Get the size of the output array from the grid height array
    dimensions=pgrheight.shape
    
    # Initialize array
    preASY=np.zeros(dimensions)
    total_weight=np.zeros(dimensions)
    
    # Calculate the asymmetry parameter per bin and sum them up to preASY
    for b in bins:    
        # Calculate the size parameter [m]
        # Set limit for empty grid cells
        size_param[b]=np.maximum(np.pi*dwet[b]/wavelength,0.001)

        # Number of particles in this bin
        n = pnaero['{}_{}'.format('NUM',b)]
        
        # Retrieve the asymmetry parameter from the look-up-table and add to total with weight n
        preASY += asymmetry(size_param[b],nr[b],ni[b],lookup)*n

        # Add up the weights
        total_weight += n
    # Finally get the mean by deviding with total weights
    ASY3D = preASY / (total_weight + 1e-12)

    # Vertical mean
    ASY = np.mean(ASY3D, axis=1)

    del ASY3D
    # beta based on [1]
    ASY2 = ASY*ASY
    ASY3 = ASY2*ASY
    term1 = ASY3*(-0.2936)
    term2 = 0.2556*ASY2
    term3 = (-0.4489)*ASY
    print("Beta is finally being calculated (all the terms are there)")
    beta_out = term1+term2
    beta_out = beta_out+term3+0.5043

    return beta_out

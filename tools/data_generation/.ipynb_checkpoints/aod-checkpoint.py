from netCDF4 import Dataset
from salsa_parameters import *
from optical_properties import cross_section
import numpy as np

def calculate_aod(nr,ni,dwet,pnaero,pgrheight,lutfile, model, wavelength):
    """
    Calculates aerosol optical depth for given volume concentration fields
    """

    # open netcdf look-up-table of optical properties
    lookup=Dataset(lutfile, 'r')

    # extinction cross section / lambda^2
    sigma=dict()
    # single scattering albedo
    omega=dict()

    # mie size parameter
    size_param=dict()
    # temporary arrays
    b_ext=dict()
    b_abs=dict()

    # calculate the extinction coefficient per bin
    for b in bins:
        
        # calculate the size parameter [m]
        # set limit for empty grid cells
        size_param[b]=np.maximum(np.pi*dwet[b]/wavelength,0.001)

        # retrieve the extinction cross section from the look-up-table
        sigma,omega=cross_section(size_param[b],nr[b],ni[b],lookup)

        # calculate extinction coefficient [1/m]
        # b_ext(bin) = sigma/lambda * N = C_ext * N
        b_ext[b]=sigma*wavelength*wavelength*pnaero['{}_{}'.format('NUM',b)]
        # calculate absorption coefficient [1/m]
        b_abs[b]=(1.0 - omega) * b_ext[b]
        # calculate absorption coefficient [1/m]
        

    # get the size of the AOD array from the grid height array
    dimensions=pgrheight.shape

    # calculate the sum of extinction coefficients  [1/m]

    aerosol_optical_depth=dict()
    aerosol_optical_depth_lev=dict()
    aerosol_absorption_optical_depth=dict()
    extinction=np.zeros(dimensions)
    absorption=np.zeros(dimensions)
    aod=np.zeros(dimensions)

    #aerosol_optical_depth=np.zeros([dimensions[0],1,dimensions[2],dimensions[3]])
    
    for b in bins:
        # calculate the extinction coefficient
        extinction+=b_ext[b]
        # calculate the absorption coefficient
        absorption+=b_abs[b]
        # calculate the aerosol optical depth
        # tau = b_ext * height of the grid box
        aod+=b_ext[b]*pgrheight



    return aod, extinction, absorption

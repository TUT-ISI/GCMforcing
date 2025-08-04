from netCDF4 import Dataset
from salsa_parameters import *
from optical_properties import cross_section
from asy_properties import asymmetry
import numpy as np


"""
by Harri Kokkola / FMI
    Atte Laakso / Aalto University

Calculates aerosol optical depth, extinction and absorption coefficients 
and average backscattering fraction for given volume concentration fields

Sources:
    [1] Taufiq Hassan, H. Moosmüller, Chul E. Chung,
        Coefficients of an analytical aerosol forcing equation determined with a Monte-Carlo radiation model,
        Journal of Quantitative Spectroscopy and Radiative Transfer,
        Volume 164, 2015, Pages 129-136, ISSN 0022-4073,
        https://doi.org/10.1016/j.jqsrt.2015.05.015.
        (https://www.sciencedirect.com/science/article/pii/S0022407315002009)
"""

def calculate_aod(nr,ni,dwet,pnaero,pgrheight,lutfile, model, wavelength):

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
    b_sca=dict()

    # variables for weighting asymmetry parameters over bins
    weighted_asymmetry_sum = np.zeros_like(pgrheight)
    total_scattering = np.zeros_like(pgrheight)

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
        # calculate scattering coefficient (b_sca = b_ext-b_abs)
        b_sca[b]=b_ext[b]-b_abs[b]
        
        # retrieve the asymmetry parameter from the look-up-table
        asym = asymmetry(size_param[b],nr[b],ni[b],lookup)

        # accumulate weighted asymmetry by scattering coefficient
        weighted_asymmetry_sum += asym * b_sca[b]

        # accumulate total scattering coefficient
        total_scattering += b_sca[b]

    # calculate average asymmetry weighted by scattering coefficient
    avg_asymmetry = np.divide(weighted_asymmetry_sum, total_scattering, out=np.zeros_like(weighted_asymmetry_sum), where=total_scattering!=0)

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

    # use scattering as weight to vertically integrate asy for beta (avg backscattering fraction)
    ASY = np.average(avg_asymmetry, weights=total_scattering, axis=1)

    # beta based on [1]
    beta = (-0.2936)*ASY**3 + 0.2556*ASY**2 + (-0.4489)*ASY + 0.5043


    print("ASY range:", np.min(ASY), np.max(ASY))
    print("β range:", np.min(beta), np.max(beta))
    print("Extinction - Absorption range:", np.min(total_scattering), np.max(total_scattering))
    print("returning beta and AOD")

    return aod, extinction, absorption, beta

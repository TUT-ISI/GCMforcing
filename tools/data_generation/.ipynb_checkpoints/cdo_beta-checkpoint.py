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


def calculate_beta(nr,ni,dwet,pnaero,pgrheight,lutfile,wavelength,selmon,tf,intime,other_model,out):
    """
    Calculates weighted average asymmetry factor for given volume concentration fields

    uses gridheight as weight in weighted mean
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

    print("3D ASY values calculated over bins, starting CDO calculations next!")
    # CDO handle the rest:
    ####################################################################
    # Utilize CDO to boost column average calculation for ASY
    os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
    cdo_path=os.getenv('CDO')
    cdo=CDO()
    cdo.setCdo(cdo_path)
    cdo.debug = True

    griddata=Dataset(pgrheight, 'r', format='NETCDF4_CLASSIC')    
    lon=griddata.variables['lon'][:]
    lat=griddata.variables['lat'][:]
    lev=griddata.variables['lev'][:]
    time=griddata.variables['time'][:]
    
    # Init some files
    ASY2 = j(tf,intime+'ASY2_temp'+other_model+'.nc')
    ASY3 = j(tf,intime+'ASY3_temp'+other_model+'.nc')
    term1 = j(tf,intime+'term1'+other_model+'.nc')
    term2 = j(tf,intime+'term2'+other_model+'.nc')
    term3 = j(tf,intime+'term3'+other_model+'.nc')
    temp1 = j(tf,intime+'totasy_temp1'+other_model+'.nc')
    temp2 = j(tf,intime+'totasy_temp2'+other_model+'.nc')
    asytemp2d = j(tf,intime+'asy_temp2d'+other_model+'.nc')
    asytemp3d = j(tf,intime+'asy_temp3d'+other_model+'.nc')
    
    print("Saving 3D ASY to netCDF")
    write_4D_grid(asytemp3d,ASY3D,lon,lat,lev,time,'beta')
    del ASY3D, preASY
    # Use CDO to get vertical average
    cdo.vertavg(
        input = ' '.join([
            selmon,
            asytemp3d,
        ]),
        output=asytemp2d
    )

    print("Calculating beta values from ASY")
    # To calculate beta using vertical average of asymmetry parameters according to polynomial in eq 5 in [1]
    
    # ASY^2
    cdo.mul( input = [asytemp2d, asytemp2d], output=ASY2 )
    
    # ASY^3
    cdo.mul( input = [ASY2, asytemp2d], output=ASY3 )
    
    # Multiply each term by its coefficient
    cdo.mulc( -0.2936, input = ASY3, output=term1 )
    cdo.mulc( 0.2556, input = ASY2, output=term2 )
    cdo.mulc( -0.4489, input = asytemp2d, output=term3 )
    
    # Add all terms together
    cdo.add( input = [term1,term2], output=temp1 )
    cdo.add( input = [temp1,term3], output=temp2 )
    cdo.addc( 0.5043, input = temp2, output=out )

    print("beta values calculated and saved to output_AOD folder")
    
    # Delete temp files
    for f in [ASY2, ASY3, temp1, temp2, term1, term2, term3, asytemp2d, asytemp3d]:
        if os.path.exists(f):
            os.remove(f)
    
    return

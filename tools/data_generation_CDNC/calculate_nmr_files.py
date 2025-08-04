"""
Written for Python3
by  Atte Laakso / UEF

Calculates the number mixing ratio using mass mixing ratio for any mmr data
 
"""

from netCDF4 import Dataset
import sys, os
import numpy as np
import netCDF4 as nc
from salsa_parameters import *
import math
from write_netcdf import *

def main():
    # SALSA path
    sp = "../SALSA_2010" # define your own
    
    # Get dry bin densities from file calculated in calculate_dry_density.py
    dry_file = os.path.join(sp,'vbs_sensitivity_base_2010.01_dens_monmean.nc')
    dry_data = nc.Dataset(dry_file, 'r', format='NETCDF4_CLASSIC')
    
    # SALSA ham
    salsa_ham = os.path.join(sp,'vbs_sensitivity_base_2010.01_ham_monmean.nc')

    tp = "./tools/AeroCom_models/" # define your own
    models = ["CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA","GFDL-AM4"]
    
    lon=dry_data.variables['lon'][:]
    lat=dry_data.variables['lat'][:]
    lev=dry_data.variables['lev'][:]
    time=dry_data.variables['time'][:]

    print(f'there are {len(bins)} bins')
    
    for i in models:
        fname_mmrham = os.path.join(tp,"aerocom3_"+i+"-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc")
        # Get PM1 for model levels
        print('Reading in PM1')
        hamdata=nc.Dataset(fname_mmrham, 'r', format='NETCDF4_CLASSIC')    
        pm1 = hamdata.variables['mmrpm1']
        mmrpm1 = np.array(pm1)

        # Since PM1 values are not binwise, we choose to weight them based on bin mean particle radius
        # as based on Athmospheric chemistry and physics by Seinfield and Pandis:
        #   "Converting from mass-based to number-based representations requires information on 
        #    particle size. Uniform mass assumptions across broad size bins can lead to 
        #    unrealistic number concentrations, particularly in ultrafine modes."
        # So weighting based on particle size (r³) leads to constraining number sizes of small aerosols
        # while still conserving total mass 

        print('Calculating nmrs')
        # Iterate over all the bins to get the weights

        # Read data from the ham file
        bindata=nc.Dataset(salsa_ham, 'r', format='NETCDF4_CLASSIC')

        # Prepare radius weights first
        radius_weights = dict()
        # Total weight sum (to normalize later)
        total_weight = np.zeros_like(mmrpm1)

        for b in bins:
            # Get the radius of the particles for the bin
            radius=bindata.variables['rdry_'+b]
            radius=np.array(radius)

            # Compute radius^3 (volume weight)
            radius = np.clip(radius, 1e-20, None)
            weight = radius**3
            
            radius_weights[b] = weight
            total_weight += weight  # add this weight to total weights

        # Calculate nmrs

        # Empty dictionary for nmr values
        number_mixing_ratio=dict()
        
        for b in bins:
            # Get the radius for this bin
            radius = radius_weights[b]**(1/3)
            
            # Calculate volume of one particle
            aerosol_volumes=4/3*math.pi*(radius**3)
            
            # Read the dry aerosol density [kg of aerosol in m³ of aerosol] of the bin
            bindensity_netc=dry_data.variables[b]
            bindensity=np.array(bindensity_netc)
            
            # Calculate mass of one aerosol particle
            aerosol_mass=aerosol_volumes*bindensity

            # Calculate the number mixing ratio
            # Avoid division by zero
            aerosol_mass[aerosol_mass < 1e-30] = np.nan
            
            # Get the weigth for bin mass and use that for nmr
            bin_mass = mmrpm1 * (radius_weights[b] / total_weight)
            nmr = bin_mass / aerosol_mass
            
            # NaN values should be zero
            nmr = np.nan_to_num(nmr, nan=0, posinf=0, neginf=0)
            nmr = np.clip(nmr, 0, None) # Limit to zero (can't be negative)

            print(f"nmr is {nmr[0,0,0,0]}")
            print(f"radius raw: {radius[0,0,0,0]} m")
            print(f"aerosol volume: {aerosol_volumes[0,0,0,0]} m³")
            print(f"bindensity: {bindensity[0,0,0,0]} kg/m³")
            print(f"aerosol_mass: {aerosol_mass[0,0,0,0]} kg")

            # Save nmr to the dicitonary under NUM-variables (for each bin)
            number_mixing_ratio['NUM_'+b]=nmr

        # Save to AeroCom files
        write_4D_grid(os.path.join(tp,"aerocom3_"+i+"-met2010_AP3-CTRL_nmr_ModelLevel_2010_monthly_T63L47.nc"),number_mixing_ratio,lon,lat,lev,time,'nmr')
        print("Ready")


if __name__ == "__main__":
    main()

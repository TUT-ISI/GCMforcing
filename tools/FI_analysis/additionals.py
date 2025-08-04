import os, sys
# define path to main folder
path_to_folder='../tools/data_generation'
# Add the module's directory to the sys.path
if path_to_folder not in sys.path:
    sys.path.append(path_to_folder)

import fnmatch
import netCDF4 as nc
import numpy as np
from write_netcdf import *
from read_tracers import read_aero_bins
from read_tracers_DRE import read_aero_bins_dp
from refractive_index import refractive_index
import math
from cdo import Cdo as CDO
from calculate_mmr import c_mmr
from salsa_parameters import bins, specs

"""
by  Atte Laakso / Aalto University

Multiple additional functions for data preparation
"""

# subroutines
def get_RH_burden(rh_file, fname_echam, vphysc_file):
    """
    Make 2D representation of RH field by calculating RH burden
    """
    # initilaize CDO environment
    os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
    cdo_path=os.getenv('CDO')
    cdo=CDO()
    
    # read temperature to a temporary netcdf file to fix indexing
    tmpfile='tmp_temp.nc'
    cdo.copy(
        input = ' '.join([
            '-sp2gp',
            '-selname,st',
            fname_echam,
        ]),
        output=tmpfile,
    )
    
    # read in gridheight data from file
    with nc.Dataset(vphysc_file, 'r', format='NETCDF4_CLASSIC') as ds:
        h = ds.variables['grheightm1']
        grid_height = np.array(h)
        
    # get relative humidity data from the file
    with nc.Dataset(rh_file, 'r', format='NETCDF4_CLASSIC') as ds:
        # naming differes between model data and salsa data
        try: 
            rh_data = ds.variables['rh']
        except:
            rh_data = ds.variables['relhum']
        relative_humidity = np.array(rh_data)   # to array

    # get temperature data from file
    with nc.Dataset(tmpfile, 'r', format='NETCDF4_CLASSIC') as ds:
        temperature_data=ds.variables['st']
        temperature = np.array(temperature_data)
    
    # remove the temporary file
    os.remove(tmpfile)
    
    # using August-Roche-Magnus approximation from Clausius–Clapeyron equation to calculate saturation vapor pressure
    temperature_C = temperature - 273.15 # convert SALSA temperature data to Celsius
    ps = 6.1094 * np.exp((17.625 * temperature_C) / (temperature_C + 243.04))
    ps = ps * 100  # hPa to Pa

    # ideal gas constant
    R = 8.3145  # (J mol-1 K-1)

    # molar mass of water
    M = 0.0180153   # (kg mol-1)

    # based on ideal gas law and definition of relative humidity, gas water content is
    gwc = (relative_humidity * ps * M) / (R * temperature)

    # to to represent this in 2D, burden must be calculated
    gwc_h = gwc*grid_height   # multipy by gridheight
    burden = np.sum(gwc_h, axis = 1)   # sum all vertical layers and get burden

    return burden


def mass_burden(mmr_file, vphysc_file):
    """
    Make 2D representation of MMR field by calculating mass burden
    """

    # read in mass mixing ratio data from file
    with nc.Dataset(mmr_file, 'r', format='NETCDF4_CLASSIC') as ds:
        ratio_data = ds.variables['mmrpm1']
        mmr = np.array(ratio_data)
        
    # read in gridheight and air density data from file
    with nc.Dataset(vphysc_file, 'r', format='NETCDF4_CLASSIC') as ds:
        h = ds.variables['grheightm1']
        grid_height = np.array(h)
        rho = ds.variables['rhoam1']
        air_density = np.array(rho)

    # multiply mmrs by corresponding heights
    mmr_z = mmr*air_density*grid_height

    # calculate the mass burden
    burden = np.sum(mmr_z, axis = 1)

    return burden


def number_burden(nmr_file, vphysc_file):
    """
    Function for calculating the number burden using files
    """

    # read in mass mixing ratio data from file
    with nc.Dataset(nmr_file, 'r', format='NETCDF4_CLASSIC') as ds:
        ratio_data = ds.variables['nmr']
        nmr = np.array(ratio_data)
        
    # read in gridheight and air density data from file
    with nc.Dataset(vphysc_file, 'r', format='NETCDF4_CLASSIC') as ds:
        h = ds.variables['grheightm1']
        grid_height = np.array(h)
        rho = ds.variables['rhoam1']
        air_density = np.array(rho)

    # multiply nmrs by corresponding heights
    nmr_z = nmr*air_density*grid_height

    # calculate the number burden
    burden = np.sum(nmr_z, axis = 1)

    return burden


def find_important(fname,indices,varname,salsa_sel):
    """
    Function for finding most important layers (highest mmr)
    """

    # importance picked maps are stored under
    imp_path = '../FI_analysis/importance_maps'
    file_name = fname.split('/')[-1]
    # check if the importance weighted map is already there
    full_path = os.path.join(imp_path, varname +'_'+ salsa_sel +'_'+ file_name)
    if os.path.exists(full_path):
        # print("Reading importance weighted data from file")
        with nc.Dataset(full_path, 'r', format='NETCDF4_CLASSIC') as ds:
            output = ds.variables[varname][:]
            output = np.array(output)
            return output
    else:
        griddata=nc.Dataset(fname, 'r', format='NETCDF4_CLASSIC')
    
        # list all available variables
        print("Available variables:", list(griddata.variables.keys()))
        
        data=griddata.variables[varname]
        # read the grid variables
        lon=griddata.variables['lon'][:]
        lat=griddata.variables['lat'][:]
        time=griddata.variables['time'][:]
    
        # initialize output list
        output=np.zeros((len(time), len(lat), len(lon)))
        # find the most effective datapoints
        for t in range(len(time)):
            for i in range(len(lat)):
                for j in range(len(lon)):
                    # get level based on index field calculated using importance_level_finder
                    ind=indices[t,i,j]
                    # update output list with most important grid value
                    output[t,i,j] = data[t,ind,i,j]
        write_3D_grid(full_path, output, lon, lat, time, varname)
        print("Saved new file to importance mapped files under parameter ",varname)
        return output


def create_refind(refrac_file, fname_tracer, fname_vphysc, real=True, model=None, mmr_file=None, mmr_comp=None):
    """
    Function for creating refractive indices files
    """

    # Volume and refractive index calculation based on given parameters
    if mmr_file and mmr_comp and model: # if mmr also included
        mratio = c_mmr(mmr_file, mmr_comp)
        pnaero, pvols = read_aero_bins_dp(fname_tracer, fname_vphysc, mratio, 'mmr')
        nr, ni = refractive_index(pvols, model)
    elif model: # If no mmr but ri
        pnaero, pvols = read_aero_bins(fname_tracer, fname_vphysc)
        nr, ni = refractive_index(pvols, model)
    else: # If no mmr and no ri
        pnaero, pvols = read_aero_bins(fname_tracer, fname_vphysc)
        nr, ni = refractive_index(pvols)

    # Read grid shape
    with nc.Dataset(fname_tracer, 'r') as griddata:
        lon = griddata.variables['lon'][:]
        lat = griddata.variables['lat'][:]
        lev = griddata.variables['lev'][:]
        time = griddata.variables['time'][:]

    # Initialize
    refind_data = np.zeros((len(time), len(lev), len(lat), len(lon)))
    weights = np.zeros_like(refind_data)

    data_dict = nr if real else ni

    for bin_name in data_dict:
        bin_ri = data_dict[bin_name]
        
        # Sum all contributions from different species
        bin_pvol = np.zeros_like(refind_data)
        for key in pvols:
            if key.endswith(f'_{bin_name}'):
                bin_pvol += pvols[key]  # sum up volume from all species in this bin

        # Weight and accumulate
        refind_data += bin_ri * bin_pvol
        weights += bin_pvol

    # Normalize
    weights[weights == 0] = np.nan  # Avoid division by zero
    refind_data /= weights

    # Save to NetCDF
    varname = 'ref_ind' if real else 'abs_coef'
    write_4D_grid(refrac_file, refind_data, lon, lat, lev, time, varname)
    return refind_data


def get_salsa_files(model, salsa_path):
    """
    Function for getting salsafile paths from the folder
    """
    # tracer file    
    fname_tracer=os.path.join(salsa_path,model+'_tracer_monmean.nc')
    # vphysc file
    fname_vphysc=os.path.join(salsa_path,model+'_vphysc_monmean.nc')
    # ham file
    fname_ham=os.path.join(salsa_path,model+'_ham_monmean.nc')
    # activ_file
    fname_activ=os.path.join(salsa_path,model+'_activ_monmean.nc')
    # echam_file
    fname_echam=os.path.join(salsa_path,model+'_echam_monmean.nc')

    return fname_tracer, fname_vphysc, fname_ham, fname_activ, fname_echam


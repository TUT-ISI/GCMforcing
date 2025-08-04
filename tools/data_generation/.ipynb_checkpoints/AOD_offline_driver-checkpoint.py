#    Harri Kokkola / FMI 
#    Atte Laakso / UEF
"""

Calculates optical properties of aerosols

"""
from netCDF4 import Dataset
from dpwet import wet_diameter
from aod import calculate_aod
from beta import calculate_beta
#from map_visuals import plot_map
import sys, os
import numpy as np
from pathlib import Path
import netCDF4 as nc
from read_tracers_CDNC import read_aero_binsdp
from read_tracers import read_aero_bins
from refractive_index import refractive_index
from write_netcdf import *
from calculate_nratio import c_nratio
from cdo import Cdo as CDO
from calculate_mmr import c_mmr
from cdnc_new import cloud_activation
from combine_levels import combine_aods, combine_abs, combine_ext
os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
cdo_path=os.getenv('CDO')

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

# Define helping subroutines

# Function for folder naming
def naming(rh,mmr,nmr,ref):
    rhname=''
    mmrname=''
    nmrname=''
    refname=''
    # Define the naming strings if spesific parameter is selected
    if rh=='rh':
        rhname='RH'
    if mmr=='mmr':
        mmrname='MMR'
    if nmr=='nmr':
        nmrname='NMR'
    if ref=='refrac':
        refname='REF_IND'
    return rhname, mmrname, nmrname, refname

# Define shorter function for creating path to file
def j(path,file):
    path=os.path.join(path,file)
    return path

# Define function for checking if similar outcomes has already been ran
def check_similiar(file_path):
    return os.path.exists(file_path)

# Define funtion for selecting month and saving data from that month to a temporary file
def cdo_copy(original_file, temporary_file, selmon):
    cdo=CDO()
    cdo.copy(
            input = ' '.join([
                selmon,
                original_file,
            ]),
            output=temporary_file
        )
    return
#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
# Function for calculating aod, abs and ext. This can be and is accessed from outside this file

def running(model,other_model,intime,light,rhsel,mmrsel,nmrsel,refsel,com_aod,avoid_save_on,paths):

    # Initial checks and definitions
    cdo=CDO()
    
    # Select month
    selmon='-seltimestep,'+intime
    
    # Output folder
    rhn, mmrn, nmrn, refn = naming(rhsel,mmrsel,nmrsel,refsel)
    output=j( j( paths[0],f"{rhn}_{mmrn}_{nmrn}_{refn}" ) , f"{other_model}_{rhn}_{mmrn}_{nmrn}_{refn}_out/" )
    
    # Ensure destination directory exists
    if not os.path.exists(output):
        os.makedirs(output)

    # SALSA model folder path
    sp=paths[1]
    # Temporary files folder
    tf=j( j( paths[2],f"temp_{rhn}_{mmrn}_{nmrn}_{refn}" ) , f"temp_{other_model}_{intime}_{rhn}_{mmrn}_{nmrn}_{refn}/" )
    # Ensure destination directory exists
    if not os.path.exists(tf):
        os.makedirs(tf)
    
    # File naming
    # Define output dimension for file naming
    if com_aod:
        dim='3D'
    else:
        dim='4D'
    # Aerosol optical depth output file
    fname2='aod_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{nmrsel}_{refsel}_for_{other_model}"+'.nc'
    # Extinction output file
    fname3='ext_2010'+intime+'_'+f"{rhsel}_{mmrsel}_{nmrsel}_{refsel}_for_{other_model}"+'.nc'
    # Absorption output file
    fname4='abs_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{nmrsel}_{refsel}_for_{other_model}"+'.nc'
    # Average backscattering coefficient file
    fname5='bac_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{nmrsel}_{refsel}_for_{other_model}"+'.nc'
    
    # Check if this iteration exists
    if avoid_save_on:
        if check_similiar(j(output,fname2)):
            if check_similiar(j(output,fname3)):
                if check_similiar(j(output,fname4)):
                    if check_similiar(j(output,fname5)):
                        # If the iteration is already completed, skip to the next one
                        print("Iteration skipped due to dublicate on file "+fname2)
                        return

    #----------------------------------------------------------------------------------------------------------
    # Read the files
    
    # Names for temporary files
    rhhamfile=j(tf,intime+'rhham_temp'+other_model+'.nc')
    mmrhamfile=j(tf,intime+'mmrham_temp'+other_model+'.nc')
    compfile=j(tf,intime+other_model+'comp_temp.nc')
    height_file=j(tf,intime+other_model+'grid_height_temp.nc')
    tracer_temp = j(tf,intime+other_model+'tracer_temp.nc')
    ham_temp = j(tf,intime+'ham_temp'+other_model+'.nc')
    
    # Tracer file 
    fname_tracer=model+'_tracer_monmean.nc'
    # Select month for tracer file
    cdo_copy( j(sp, fname_tracer) , tracer_temp, selmon)
    tracerfile=tracer_temp # Reallocate the selected data back to original variable

    # vphysc file
    fname_vphysc=model+'_vphysc_monmean.nc'
    # Select month for vphysc file and copy that data to a temporary file
    cdo_copy( j(sp, fname_vphysc), height_file, selmon)
    vphyscfile=height_file
    
    # ham_file
    fname_ham = model+'_ham_monmean.nc'
    # Select month for the data and copy data to file
    cdo_copy( j(sp,fname_ham) , ham_temp, selmon)
    fname_ham=ham_temp
    
    # ham_file for relative humidity (RH)
    if rhsel=='rh':   # select the model data only if RH field is to be changed
        fname_rhham='aerocom3_'+other_model+'-met2010_AP3-CTRL_rh_ModelLevel_2010_monthly_T63L47.nc'
        fname_rhham=j(paths[3],fname_rhham) # Find correct path to file
        # Select month an copy data to a temporary file
        cdo_copy(fname_rhham, rhhamfile, selmon)
        fname_rhham = rhhamfile
    else:
        # Use SALSA ham_file (initialized above)
        fname_rhham = ham_temp
    
    # File for mass mixing ratio (MMR)
    fname_mmrham='aerocom3_'+other_model+'-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc'
    fname_mmrham=j(paths[3],fname_mmrham) # Find correct path to file
    # Select time for MMR data
    cdo_copy(fname_mmrham, mmrhamfile, selmon)
    fname_mmrham = mmrhamfile
    
    # File for MMR comparison
    fname_mmrcomp='aerocom3_ECHAM6.3-SALSA2.0-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc'
    fname_mmrcomp=j(paths[3],fname_mmrcomp) # Find correct path to file
    # Select month for MMR comparison data
    cdo_copy(fname_mmrcomp, compfile, selmon)
    fname_mmrcomp = compfile

    # Check if files exist 
    if not os.path.exists(tracerfile):
        print("Tracer file not found.")        
        exit()

    if not os.path.exists(vphyscfile):
        print("vphysc file not found.")
        exit()
        
    #----------------------------------------------------------------------------------------------------------
    # Handle the data from files
    
    # Wavelength in nm
    lmbd_str=light
    wavelength=int(lmbd_str)*1e-9
    
    relhum_fixed=-999

    # String for a fixed rh (e.g. rh80 for 80%)
    rh_str='rh'+str(int(relhum_fixed*100))

    # If RH=0, hydration is not calculated
    if relhum_fixed==0:
        rh_str='dry'

    # If ambient model RH is used
    if relhum_fixed<0:
        rh_str=''
    
    # Load grid data from file
    griddata=nc.Dataset(vphyscfile, 'r', format='NETCDF4_CLASSIC')    
    lon=griddata.variables['lon'][:]
    lat=griddata.variables['lat'][:]
    lev=griddata.variables['lev'][:]
    time=griddata.variables['time'][:]

    # Define if outcome for aod and abs is 3D or 4D
    if (com_aod):
        griddim=3
    else:
        griddim=4

    # Calculate the mass mixing ratios' ratio between other_model and SALSA
    ratio=c_mmr(fname_mmrham,fname_mmrcomp)

    # dp selection for read_aero_bins function
    if "salsa" in other_model.lower():
        dp=1
    else:
        dp=0
        
    # Get number and volume concentrations of individual species
    print('Reading in number and volume concentrations')
    # Use recalculated bins if mmr or nmr fields are changed
    if mmrsel=='mmr' or nmrsel=='nmr':
        znaero,zvols = read_aero_binsdp(tracerfile,vphyscfile,fname_mmrham,fname_ham,lon,lat,ratio,mmrsel,nmrsel)
        # Clear no longer used variables
        del ratio
    else:
        znaero,zvols = read_aero_bins(tracerfile,vphyscfile,lon,lat)  # Otherwise use SALSA bins

    # Get grid cell height for calculating AOD
    print('Reading in grid box height from vphysc_file.nc')
    vphyscdata = Dataset(vphyscfile, 'r')
    h = vphyscdata.variables['grheightm1']
    zgrheight = np.array(h)

    print(f"gridheight shape is {zgrheight.shape}")
    
    # Get relative humidity for model levels
    print('Reading in relative humidity from rhham_file.nc')
    if rhsel=='rh':
        rhhamdata = nc.Dataset(fname_rhham, 'r', format='NETCDF4_CLASSIC')   # Use RH data from the other model    
        relhum = rhhamdata.variables['rh']
    else:
        rhhamdata = nc.Dataset(fname_rhham, 'r', format='NETCDF4_CLASSIC')    # Use RH data from SALSA
        relhum = rhhamdata.variables['relhum']
    zrh = np.array(relhum)   # Relative humidity values

    print(f"relative humidity shape is {zrh.shape}")
    
    # Clear unused variables
    try:
        rhhamdata.close()
    except IOError as e:
        pass
    del rhhamdata, relhum   
    
    
    # Calculate the wet diameter of particles in each size bin
    print('Calculating the wet diameter')
    dwet,ddry = wet_diameter(znaero,zvols,zrh,model,relhum_fixed)
    
    # Calculate the refractive index based on volume concentrations
    print('Calculating refractive indices')
    # If refractive indices are changed to the ones of the model, read them from another file
    if refsel=='refrac':
        nr,ni=refractive_index(zvols, other_model)
    # Else read indices from SALSA's parameters
    else:
        nr,ni=refractive_index(zvols)
    
    print('Calculating AOD')
    lut_file='lut_optical_properties.nc'
    aod,extinction,absorption=calculate_aod(nr,ni,dwet,znaero,zgrheight,lut_file, model, wavelength)
    
    print('Calculating average backscattering coefficient')
    # calculate_beta(nr,ni,dwet,znaero,zgrheight,lut_file, wavelength, selmon, tf, intime, other_model, j(output,fname5))
    
    # other version
    beta = calculate_beta(nr,ni,dwet,znaero,zgrheight,lut_file, wavelength)
    write_3D_grid(j(output,fname5),beta,lon,lat,time,'beta')
    
    # Remove the temporary files and directory
    if rhsel=='rh':
        os.remove(rhhamfile)
    os.remove(mmrhamfile)
    os.remove(compfile)
    os.remove(height_file)
    os.remove(tracer_temp)
    os.remove(ham_temp)
    os.rmdir(tf)
    #----------------------------------------------------------------------------------------------------------
    # Save data to files
    
    print('Saving data to NetCDF files')


    if(griddim == 2):
        
        write_2D_grid('dwet.nc',dwet,lon,lat)
    
        write_2D_grid('ddry.nc',ddry,lon,lat)
        
        write_2D_grid('nr.nc',nr,lon,lat)
        
        write_2D_grid('ni.nc',ni,lon,lat)
        
        write_2D_grid('rh.nc',rh,lon,lat)
        
        write_2D_grid(j(output,fname2),aerosol_optical_depth,lon,lat,'AOD')

    
    elif(griddim == 3):
        # Combine all output levels to 2D if option selected (if griddim=3)
        if com_aod:
            aod = combine_aods(aod,lat,lon,time)
            absorption = combine_abs(absorption,lat,lon,time)
            extinction = combine_ext(extinction,lat,lon,time)

        # Write netCDFs
        write_3D_grid(j(output,fname2),aod,lon,lat,time,'AOD')
        write_3D_grid(j(output,fname3),extinction,lon,lat,time,'ec'+lmbd_str+rh_str+'aer')
        write_3D_grid(j(output,fname4),absorption,lon,lat,time,'ac'+lmbd_str+rh_str+'aer')
    
    elif(griddim == 4):
        
        write_4D_grid(j(output,fname2),aod,lon,lat,lev,time,'AOD')
        write_4D_grid(j(output,fname3),extinction,lon,lat,lev,time,'ec'+lmbd_str+rh_str+'aer')
        write_4D_grid(j(output,fname4),absorption,lon,lat,lev,time,'ac'+lmbd_str+rh_str+'aer')

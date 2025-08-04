from netCDF4 import Dataset
from dpwet import wet_diameter
from aod import calculate_aod
#from map_visuals import plot_map
import sys, os
import numpy as np
from pathlib import Path
import netCDF4 as nc
from read_tracers_DRE import read_aero_bins_dp
from read_tracers import read_aero_bins
from refractive_index import refractive_index
from write_netcdf import *
from cdo import Cdo as CDO
from calculate_mmr import c_mmr
from combine_levels import combine_aods, combine_abs
import time as t
os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
cdo_path=os.getenv('CDO')

"""
by  Harri Kokkola / FMI 
    Atte Laakso / Aalto University

Calculates optical properties of aerosols and aerosol direct radiative effect

Sources:
[1] Variability of Aerosol Optical Properties at Four North American Surface Monitoring Sites (2002)
    David J. Delene and John A. Ogren
    Page(s): 1135-1150
    link to the article: https://journals.ametsoc.org/view/journals/atsc/59/6/1520-0469_2002_059_1135_voaopa_2.0.co_2.xml?tab_body=pdf
    doi: https://doi.org/10.1175/1520-0469(2002)059<1135:VOAOPA>2.0.CO;2
"""

# Define helping subroutines

# Function for folder naming
def naming(rh,mmr,ref,clt):
    rhname=''
    mmrname=''
    refname=''
    cltname=''
    # Define the naming strings if spesific parameter is selected
    if rh=='rh':
        rhname='RH'
    if mmr=='mmr':
        mmrname='MMR'
    if ref=='refrac':
        refname='REF_IND'
    if clt=='clt':
        cltname='CLT'
    return rhname, mmrname, refname, cltname

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
# Function for calculating aod, abs, beta and ADRE. This can be and is accessed from outside this file

def running(model,other_model,intime,light,rhsel,mmrsel,refsel,cltsel,avoid_save_on,paths):
    # Keep track of time
    start_time = t.time()

    # Initial checks and definitions
    cdo=CDO()
    
    # Select month
    selmon='-seltimestep,'+intime
    
    # Output folder
    rhn, mmrn, refn, cltn = naming(rhsel,mmrsel,refsel,cltsel)
    output=j( j( paths[0],f"{rhn}_{mmrn}_{refn}_{cltn}" ) , f"{other_model}_{rhn}_{mmrn}_{refn}_{cltn}_out/" )
    
    # Ensure destination directory exists
    if not os.path.exists(output):
        os.makedirs(output)

    # SALSA model folder path
    sp=paths[1]
    # Temporary files folder
    tf=j( j( paths[2],f"temp_{rhn}_{mmrn}_{refn}_{cltn}" ) , f"temp_{other_model}_{intime}_{rhn}_{mmrn}_{refn}_{cltn}/" )
    # Ensure destination directory exists
    if not os.path.exists(tf):
        os.makedirs(tf)
    
    # File naming
    # Define output dimension for file naming
    dim='3D'
    # Aerosol optical depth output file
    fname2='aod_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{refsel}_{cltsel}_for_{other_model}"+'.nc'
    # Aerosol single scattering albedo output file
    fname3='ssa_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{refsel}_{cltsel}_for_{other_model}"+'.nc'
    # Absorption output file
    fname4='abs_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{refsel}_{cltsel}_for_{other_model}"+'.nc'
    # Average backscattering coefficient file
    fname5='bac_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{refsel}_{cltsel}_for_{other_model}"+'.nc'
    # Aerosol direct radiative effect file
    fname6='ADRE_2010'+intime+'_'+dim+f"{rhsel}_{mmrsel}_{refsel}_{cltsel}_for_{other_model}"+'.nc'

    # Check if this iteration exists
    if avoid_save_on:
        if check_similiar(j(output,fname2)):
            if check_similiar(j(output,fname3)):
                if check_similiar(j(output,fname4)):
                    if check_similiar(j(output,fname5)):
                        if check_similiar(j(output,fname6)):
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
    echam_temp = j(tf,intime+'echam_temp'+other_model+'.nc')
    clt_temp = j(tf,intime+'clt_temp'+other_model+'.nc')

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

    # echam file
    fname_echam = model+'_echam_monmean.nc'
    # Select month for the data and copy data to file
    cdo_copy( j(sp,fname_echam) , echam_temp, selmon)
    fname_echam=echam_temp

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

    # File for total cloud fraction
    fname_clt='aerocom3_'+other_model+'-met2010_AP3-CTRL_clt_ModelLevel_2010_monthly_T63L47.nc'
    fname_clt=j(paths[3],fname_clt) # Find correct path to file
    # Select time for cloud fraction data
    cdo_copy(fname_clt, clt_temp, selmon)
    fname_clt = clt_temp

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

    # Calculate the mass mixing ratios' ratio between other_model and SALSA
    ratio=c_mmr(fname_mmrham,fname_mmrcomp)
        
    # Get number and volume concentrations of individual species
    print('Reading in number and volume concentrations')
    # Use recalculated bins if mmr fields are changed
    if mmrsel=='mmr':
        znaero,zvols = read_aero_bins_dp(tracerfile,vphyscfile,ratio,mmrsel)
        # Clear no longer used variables
        del ratio
    else:
        znaero,zvols = read_aero_bins(tracerfile,vphyscfile)  # Otherwise use SALSA bins

    # Get grid cell height for calculating AOD
    print('Reading in grid box height from vphysc_file.nc')
    vphyscdata = Dataset(vphyscfile, 'r')
    h = vphyscdata.variables['grheightm1']
    zgrheight = np.array(h)
    
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
    
    #--------------------------------------------------------------------------------------------------
    # Central calculations

    print('Calculating AOD')
    lut_file='lut_optical_properties.nc'
    AOD,extinction,absorption,backscatter = calculate_aod(nr,ni,dwet,znaero,zgrheight,lut_file, model, wavelength)
    
    # Load grid data from file
    griddata=nc.Dataset(vphyscfile, 'r', format='NETCDF4_CLASSIC')    
    lon=griddata.variables['lon'][:]
    lat=griddata.variables['lat'][:]
    lev=griddata.variables['lev'][:]
    time=griddata.variables['time'][:]

    print('Reading in CLT and albedo data')
    with nc.Dataset(fname_echam, 'r', format='NETCDF4_CLASSIC') as sf:
        alb_var = sf.variables['albedo']
        albedo = np.array(alb_var)

        clt_var = sf.variables['aclcov']
        clouds = np.array(clt_var)
    # If clouds are selected, then read the cloud fraction data from other model
    if cltsel=='clt':
        with nc.Dataset(fname_clt, 'r', format='NETCDF4_CLASSIC') as cf:
            clt_var = cf.variables['clt']
            clouds = np.array(clt_var)
            # GISS and 
            if "GISS" in other_model:
                clouds = clouds / 10000 # GISS cloud fraction is in units of 100th of a %
            elif "GFDL" in other_model:
                clouds = clouds / 100   #GFDL has cloud fractions as % so devide by 100
            if clouds.ndim == 4:
                print("Detected 4D cloud fraction data. Estimating cover over all levels.")
                # Use maximum overlap approximation for calculating total cover
                clouds = np.max(clouds, axis=1)

    # Cloud fraction must be between 0 and 1
    clouds = np.clip(clouds, 0, 1)
    print(f"Average cloud cover is {np.mean(clouds)}")

    del alb_var, clt_var

    print("Combining levels")
    # Combine all output levels to 3D
    AOD = combine_aods(AOD,lat,lon,time)
    # Multiply 4D values of absorption with grid heights
    # to get the absorption optical depth values for SSA (omega) calculation (in c_DRE)
    absorption = combine_abs(absorption*zgrheight,lat,lon,time)

    print("Calculating ADRE")
    
    # Method for calculating the incoming direct radiative effect using 
    # other values and based on eq (6) in [1]
    def c_DRE(tau,tau_abs,beta, Ac, Rs):
        # Initialize constants
        D = 0.5           # the fractional day length
        S_0 = 1370        # the solar constant (W/m-2)
        T_at = 0.76       # the atmospheric transmission
        # By definition
        omega_0 = 1-tau_abs/tau
        # eq (6) in [1]
        DRE = -D*S_0*(T_at**2)*(1-Ac)*tau * ( beta*omega_0*(1-Rs)**2 - 2*(1-omega_0)*Rs )
        return DRE, omega_0

    # Calculate ADRE
    ADRE, omega = c_DRE(AOD,absorption,backscatter,clouds,albedo)

    # Remove the temporary files and directory
    if rhsel=='rh':
        os.remove(rhhamfile)
    os.remove(mmrhamfile)
    os.remove(compfile)
    os.remove(height_file)
    os.remove(tracer_temp)
    os.remove(ham_temp)
    os.remove(clt_temp)
    if os.path.isdir(tf) and not os.listdir(tf):
        os.rmdir(tf)

    #----------------------------------------------------------------------------------------------------------
    # Save data to files
    
    print('Saving data to NetCDF files')
    # Write netCDFs
    write_3D_grid(j(output,fname2),AOD,lon,lat,time,'AOD')
    write_3D_grid(j(output,fname3),omega,lon,lat,time,'SSA')
    write_3D_grid(j(output,fname4),absorption,lon,lat,time,'ac'+lmbd_str+rh_str+'aer')
    write_3D_grid(j(output,fname5),backscatter,lon,lat,time,'beta')
    write_3D_grid(j(output,fname6),ADRE,lon,lat,time,'ADRE')

    total_time = t.time()-start_time
    print(f"Iteration with {rhsel}, {mmrsel}, {refsel}, {cltsel} completed in {total_time}")


# Written for Python v2.7.6
# by Harri Kokkola / FMI 
#  & Atte Laakso / UEF
"""

Calculates optical properties 
 
"""
from netCDF4 import Dataset
from dpwet import wet_diameter
#from map_visuals import plot_map
import sys, os, shutil
import numpy as np
import math
from pathlib import Path
import netCDF4 as nc
from read_tracers_CDNC import read_aero_binsdp
from read_tracers import read_aero_bins
from write_netcdf import *
from calculate_nratio import c_nratio
from cdo import Cdo as CDO
from calculate_mmr import c_mmr
from cdnc_new_dp import cloud_activation_dp
os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
cdo_path=os.getenv('CDO')

#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
# define helping subroutines
    
# function for folder naming
def naming(mmr,nmr,updraft,clt,clw):
    mmrname=''
    nmrname=''
    updraftname=''
    cltname=''
    clwname=''
    # define the naming strings if spesific parameter is selected
    if mmr=='mmr':
        mmrname='MMR'
    if nmr=='nmr':
        nmrname='NMR'
    if updraft=='updraft':
        updraftname='UPD'
    if clt=='clt':
        cltname='CLT'
    if clw=='clw':
        clwname='CLW'
    return mmrname, nmrname, updraftname, cltname, clwname

#define shorter function for creating path to file
def j(path,file):
    path=os.path.join(path,file)
    return path

#define function for checking if similar outcomes has already been ran
def check_similiar(file_path):
    return os.path.exists(file_path)
    
# define funtion for selecting month and saving data from that month to a temporary file
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
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

# code for calculating CDNC and energy balance with aerosol cloud radiative effect

# define the function for operating this code outside this file
def running(model,other_model,intime,light,mmrsel,nmrsel,updraftsel,cltsel,clwsel,avoid_save_on,calculate_CRE):
    
    #--------------------------------------------------------------------------------------------------------------------------------
    # set general commands and definitions
    
    cdo=CDO()    

    # command for selecting a month
    selmon='-seltimestep,'+intime
    
    # output folder
    mmrn, nmrn, updn, cltn, clwn = naming(mmrsel,nmrsel,updraftsel,cltsel,clwsel)
    output='CDNC_output/'+f"{mmrn}_{nmrn}_{updn}_{cltn}_{clwn}/"+f"{other_model}_{mmrn}_{nmrn}_{updn}_{cltn}_{clwn}_out/"
    
    # ensure destination directory exists
    if not os.path.exists(output):
        os.makedirs(output)
        
    # SALSA model folder path
    sp='/scratch/project_2010692/SALSA_2010'
    # temporary files folder
    tf='temporary/'+f"temp_{mmrn}_{nmrn}_{updn}_{cltn}_{clwn}/"+f"temp_{other_model}_{intime}_{mmrn}_{nmrn}_{updn}_{cltn}_{clwn}/"
    
    # ensure destination directory exists
    if not os.path.exists(tf):
        os.makedirs(tf)
    
    # file naming
    
    dim='4D'   # define output dimension for file naming
    
    # cloud activation output file
    fname5='cdnc_2010'+intime+'_'+dim+f"{mmrn}_{nmrn}_{updn}_{cltn}_{clwn}_for_{other_model}"+'.nc'
    fname6='CRE_2010'+intime+'_'+dim+f"{mmrn}_{nmrn}_{updn}_{cltn}_{clwn}_for_{other_model}"+'.nc'

    # check if this iteration exists
    if avoid_save_on:
        if calculate_CRE:
            if check_similiar(output+fname6):
                # if the iteration is allready completed, skip to the next one
                print("Iteration skipped due to dublicate on file "+fname6)
                return
        else:
            if check_similiar(output+fname5):
                # if the iteration is allready completed, skip to the next one
                print("Iteration skipped due to dublicate on file "+fname5)
                return
    #--------------------------------------------------------------------------------------------------------------------------------
    # set files and select month
    
    # names for temporary files
    mmrhamfile=j(tf,intime+'mmrham_temp'+other_model+'.nc')
    tmpfile=j(tf,'tmp_temp'+intime+other_model+'.nc')
    compfile=j(tf,intime+other_model+'comp_temp.nc')
    tracer_temp = j(tf,intime+other_model+'tracer_temp.nc')
    vphysc_temp = j(tf,intime+'vphysc_temp'+other_model+'.nc')
    activ_temp = j(tf,intime+'activ_temp'+other_model+'.nc')
    echam_temp = j(tf,intime+'echam_temp'+other_model+'.nc')
    ham_temp = j(tf,intime+'ham_temp'+other_model+'.nc')
    
    # tracer file 
    fname_tracer = model+'_tracer_monmean.nc'
        # select month for tracer file
    cdo_copy( j(sp, fname_tracer) , tracer_temp, selmon)
    tracerfile = tracer_temp
    
    # vphysc file
    fname_vphysc = model+'_vphysc_monmean.nc'
        # select month for vphysc file and copy that data to a temporary file
    cdo_copy( j(sp, fname_vphysc), vphysc_temp, selmon)
    vphyscfile = vphysc_temp

    # ham_file
    fname_ham = model+'_ham_monmean.nc'
    # select month for the data and copy data to file
    cdo_copy( j(sp,fname_ham) , ham_temp, selmon)
    fname_ham=ham_temp
    
    # other files' month selections and definitions
    
    # for mass mixing ratio (MMR)
    fname_mmrham = 'aerocom3_'+other_model+'-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc'
        # select month for MMR data
    cdo_copy(fname_mmrham, mmrhamfile, selmon)
    fname_mmrham = mmrhamfile
    
    # for MMR comparison
    fname_mmrcomp = 'aerocom3_ECHAM6.3-SALSA2.0-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc'
        # select month for MMR comparison data
    cdo_copy(fname_mmrcomp, compfile, selmon)
    fname_mmrcomp = compfile

    # activ_file
    fname_activ = model+'_activ_monmean.nc'
        # select month to activ file
    cdo_copy(j(sp,fname_activ), activ_temp, selmon)
    fname_activ = activ_temp

    # echam_file
    fname_echam=model+'_echam_monmean.nc'
        # select month to echam file
    cdo_copy(j(sp,fname_echam), echam_temp, selmon)
    fname_echam = echam_temp    

    # read temperature to a temporary netcdf file
    cdo.copy(
        input = ' '.join([
            '-sp2gp',
            '-selname,st',
            fname_echam,
        ]),
        output=tmpfile,
    )

    #--------------------------------------------------------------------------------------------------------------------------------
    # do the calculations for CDNC
    
    # wavelength in nm
    lmbd_str=light
    wavelength=int(lmbd_str)*1e-9
            
    # check if files exist 
    if not os.path.exists(tracerfile):
        print("Tracer file not found.")        
        exit()

    if not os.path.exists(vphyscfile):
        print("vphysc file not found.")
        exit()

    # load grid data from file
    griddata=nc.Dataset(vphyscfile, 'r', format='NETCDF4_CLASSIC')    
    lon=griddata.variables['lon'][:]
    lat=griddata.variables['lat'][:]
    lev=griddata.variables['lev'][:]
    time=griddata.variables['time'][:]

    # calculate the ratio between mass mixing ratios of other_model and SALSA
    ratio=c_mmr(fname_mmrham,fname_mmrcomp)
       
    # get number and volume concentrations of individual species
    print('Reading in number and volume concentrations')
    
    # use recalculated bins if mmr or nmr fields are changed (bins will be scaled accordingly)
    if mmrsel=='mmr' or nmrsel=='nmr':
        znaero,zvols = read_aero_binsdp(tracerfile , vphyscfile , fname_mmrham , fname_ham , lon , lat , ratio , mmrsel , nmrsel)
    else:
        znaero,zvols = read_aero_bins(tracerfile ,  vphyscfile , lon , lat)

    # calculate full-level pressure (pfull) from half-level values (phalf)
    vphyscdata = Dataset(vphyscfile, 'r')
    phalf = vphyscdata.variables['aphm1']
    pfull = (phalf[:,0:47,:,:] + phalf[:,1:48,:,:])/2.
    zapm1 = np.array(pfull)
    # clear used variables
    del phalf, pfull
    
    # get updraft velocities for model levels
    print('Reading in updraft velocities')
    if updraftsel=='updraft':
        fname_updraft='aerocom3_'+other_model+'-met2010_AP3-CTRL_w_ModelLevel_2010_monthly_T63L47.nc'   # define file name
        
        # select month for updraft and copy data to temporary file
        updraft_temp = j(tf,intime+'updraft_temp'+other_model+'.nc')
        cdo_copy(fname_updraft,updraft_temp,selmon)
        
        activdata=nc.Dataset(updraft_temp, 'r', format='NETCDF4_CLASSIC')   # get data from file
        vervel = activdata.variables['W']   # read updraft from variables
        zw = np.array(vervel)
        os.remove(updraft_temp)   # remove updraft temporary file
    else: 
        activdata=nc.Dataset(fname_activ, 'r', format='NETCDF4_CLASSIC')   # get updraft data from SALSA files
        vervel = activdata.variables['W']
        zw = np.array(vervel)
    
    print('Reading maximum supersaturation')
    swat_max_strat = activdata.variables['SWAT_MAX_STRAT']
    psmax=swat_max_strat
    
    # delete unused variables
    try:
        activdata.close()
    except IOError as e:
        pass
    del activdata, vervel, swat_max_strat
    
    # get specific humidity for model levels
    echamdata=nc.Dataset(fname_echam, 'r', format='NETCDF_CLASSIC')
    specific_humidity=echamdata.variables['q']
    zqm1 = np.array(specific_humidity)
    
    # get temperature for model levels
    print("Reading in temperature for model levels.")
    temperaturedata=nc.Dataset(tmpfile, 'r', format='NETCDF_CLASSIC')
    temperature=temperaturedata.variables['st']
    ztm1 = np.array(temperature)
    
    # remove unused variables
    try:
        temperaturedata.close()
    except IOError as e:
        pass
    del specific_humidity, temperaturedata, temperature
    
    # remove the temporary files and directory
    if updraftsel=='updraft':
        os.remove(updraft_temp)
    os.remove(tmpfile)
    os.remove(mmrhamfile)
    os.remove(compfile)
    os.remove(tracer_temp)
    os.remove(vphysc_temp)
    os.remove(activ_temp)
    os.remove(echam_temp)
    if not calculate_CRE:
        os.rmdir(tf)

    # number of updrafts
    nw=1

    print('Calculating number of activated droplets')
    cdnc,zsmax=cloud_activation_dp(znaero, zvols, ztm1, zapm1, zqm1, zw, nw, psmax, 0, ratio)
    
    print('Saving CDNC data to NetCDF files')
    # store data to NetCDF file
    write_4D_grid(output+fname5,cdnc,lon,lat,lev,time,'CDNC')

    #--------------------------------------------------------------------------------------------------------------------------------
    # do the calculations for CRE if needed
    
    # if cloud radiative effect is wanted, calcualte the energy balance with that effect based on CDNC
    if calculate_CRE:
        print('Calculating aerosol cloud radiative effect')

        # equations based on article 'A simple model of global aerosol indirect effects'
        # link to the article:  https://doi.org/10.1002/jgrd.50567
       
        print("Reading in cloud fraction data for model levels.")
        # determine cloud fraction f_c
        if cltsel=='clt':
            # read cloud fraction data from model data if selected
            cltfile='aerocom3_'+other_model+'-met2010_AP3-CTRL_clt_ModelLevel_2010_monthly_T63L47.nc'
            
            # select month for the data
            clt_temp = j(tf,intime+'clt_temp'+other_model+'.nc')
            cdo_copy(cltfile, clt_temp, selmon)

            # use month-selected data to get cloud fraction values
            with nc.Dataset(clt_temp, 'r', format='NETCDF4_CLASSIC') as ds:
                fc_data=ds.variables['clt'][:]
                f_c=np.array(fc_data)
        else:
            # read cloud fraction data from SALSA data
            fc_data=echamdata.variables['aclcac'][:]
            f_c=np.array(fc_data)            
            
        # get liquid water content
        if clwsel=='clw':
            # read cloud fraction data from model data if selected
            clwfile='aerocom3_'+other_model+'-met2010_AP3-CTRL_clw_ModelLevel_2010_monthly_T63L47.nc'

            # select month for the data
            clw_temp = j(tf,intime+'clw_temp'+other_model+'.nc')
            cdo_copy(clwfile, clw_temp, selmon)
            
            # use month-selected data to get cloud liquid water content
            with nc.Dataset(clw_temp, 'r', format='NETCDF4_CLASSIC') as ds:
                clw_data = ds.variables['clw']
                clw = np.array(clw_data)
        else:
            # read cloud fraction data from SALSA data
            clw_data = echamdata.variables['xl']
            clw = np.array(clw_data)

        print(f"Cloud liquid water content dimensions are {clw.shape}")

        # read air density
        rho_a = vphyscdata.variables['rhoam1']
        rho_a = np.array(rho_a)
        
        # function for finding right levels for cloud coverage by selecting the lowest level where clw is greater than 0,01 g/m^3
        def find_cloud_cover_levels(clw, rho_a, f_c):
            # transform clw from kg/kg to g/m^3 by multiplying with air density and by multiplying with 1000
            lwc = clw * rho_a * 1000
            
            # define shape of the data
            month, levels, lat, lon = lwc.shape
            
            # initialize array for indices
            indices = np.full((month,lat,lon) , -1)

            # iterate trough each gridpoint
            for m in range(month):
                for i in range(lat):
                    for j in range(lon):
                        # find the all the indices of lwc values greater than 0.01 g/m^3
                        valid_indices = np.where(lwc[m,:,i,j] > 0.01)[0]
                        
                        # add only if there was a value over 0,01 g/m^3
                        if valid_indices.size != 0:
                            indices[m,i,j] = valid_indices[-1]   # highest index = lowest level
             
            # find the cloud cover levels based on lowest levels where the cloud liquid water content is over 0,01 g/m^3
            cloud_cover = np.zeros((month,lat,lon))   # initialize array
            
            # iterate over grid points
            for m in range(month):
                for i in range(lat):
                    for j in range(lon):
                        index = indices[m,i,j]   # get the corresponding index
                        if index != -1:
                            cloud_cover[m,i,j] = f_c[m,index,i,j]   # get right cloud cover value
                            # if there was no value over 0,01 g/m^3, cloud cover remains zero
                            
            return cloud_cover

        # form cloud cover
        # cloud cover values are read from the lowest level, where the lwc value exeeds 0,01 g/m^3
        # if the condition is not met, cloud cover is zero
        cloud_cover = find_cloud_cover_levels(clw, rho_a, f_c)

        print(f'Mean cloud cover is {np.mean(cloud_cover)}')
        
        # get grid box height
        grid_height = vphyscdata.variables['grheightm1']
        grid_height = np.array(grid_height)

        # calculate properties needed for cloud albedo
        print("Calculating properties for cloud albedo.")

        def calculate_effective_radius(lwc, nd, rho_w=1000):
            
            re = (3 * lwc * rho_a / (4 * np.pi * rho_w * nd)) ** (1.0 / 3.0)  # lwc is per unit mass so multiplication with density required
            
            # round small values (and too small nd value spots that are now infinity) to zero
            re = np.nan_to_num(re, nan=0, posinf=0, neginf=0)
                        
            return re
            
        def calculate_liquid_water_path(lwc_profile, dz):
            # (here also profile is per unit mass so multiplication with density required as well as with 1000 as grams should be used instead of kilograms)
            return np.sum(lwc_profile * rho_a * 1000 * dz, axis=1)
        
        def calculate_cloud_optical_thickness(lwp, re, rho_w=1000):
            tau = (3 / 2) * (lwp / (re * rho_w  * 1000))   # change density to g/m^3
            
            # replace NaN values (also: infinity = no effective radius = no optical thickness)
            tau = np.nan_to_num(tau, nan=0, posinf=0, neginf=0)
            
            return tau
        
        def calculate_cloud_albedo(tau):
            return tau / (tau + 7.7)
        
        def calculate_cloud_albedo_from_profiles(lwc_profile, nd_profile, dz):
            # use functions defined above and calculate cloud albedo
            
            re_profile = calculate_effective_radius(lwc_profile, nd_profile)   # in meters
            
            lwp = calculate_liquid_water_path(lwc_profile, dz)  # in g/m^2
            
            print("Mean liquid water path is (g/^2)",np.mean(np.mean(lwp)))
            
            re_mean = np.mean(re_profile[ re_profile > 0 ])   # in meters (uses only positive (existing) values)
            print("Mean effective radius is ",re_mean)
            
            tau = calculate_cloud_optical_thickness(lwp, re_mean)   # unitless
            
            albedo = calculate_cloud_albedo(tau)   # unitless
            
            return albedo

        # set profiles and grid heights
        lwc_profile = clw    # kg/kg
        dz = grid_height     # grid spacings in meters
        nd_profile = cdnc    # 1/m^3 
        
        # calculate cloud albedo
        cloud_albedo = calculate_cloud_albedo_from_profiles(lwc_profile, nd_profile, dz)
        
        print("Global mean cloud albedo is ",np.mean(cloud_albedo))
        
        # calculate the energy balance
        print("Calculating energy balance with cloud radiative radiative effect")

        S_0 = 1367   # downward solar at the top of the atmosphere (W m-2)

        as_data = echamdata.variables['albedo'][:]   # read surface albedo from salsa data
        surface_albedo = np.array(as_data)

        print(f"Mean surface albedo is {np.mean(surface_albedo)}")
        
        # calculate energy balance with aerosols
        E = 0.25 * S_0 * (
            (1 - cloud_cover) * (1 - surface_albedo) +
            ((cloud_cover * (1 - cloud_albedo) * (1 - surface_albedo)) / (1 - cloud_albedo * surface_albedo))
            )

        print(f"E dimensions are {E.shape}")

        # delete temporary files and folders if they still exist     
        if os.path.isdir(tf):
            shutil.rmtree(tf)

        print("Saving data to files")
        # write data to NetCDF-file
        if os.path.exists(output+fname6):  
            os.remove(output+fname6)   # if the file exists already, delete it   
        write_3D_grid(output+fname6,E,lon,lat,time,'E')

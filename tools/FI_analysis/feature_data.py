import numpy as np
import netCDF4 as nc
import os, sys
import math
from additionals import *

"""
by Atte Laakso / Aalto University

Create and return key atmospheric parameter grids for a given model configuration.

Parameters:
- parameters: list or dict indicating which parameters are to be read from model or fallback data
- model: string identifier for the model to use
Returns:
- Dictionary with keys ['rh', 'mmr', 'refrac', 'abs_coef', 'clt']
"""

# helping function
def read_netcdf_variable(path, variable, level):
    with nc.Dataset(path, 'r', format='NETCDF4_CLASSIC') as ds:
        data = np.array(ds.variables[variable])
        if data.ndim == 4:
            return data[:, level, :, :]
        return data


# actual function for returning dictionary with feature data 
def create_combos(parameters, model, path_to_folder, fname_ham, fname_echam, fname_vphysc, fname_tracer):
    
    print('Initializing datasets for training and test data for model '+model+'.')

    # library to store the data
    root_data={}
    
    # initialize all the variables to be stated
    rh_data=[]
    mmr_data=[]
    refractive_ind_data=[]
    abs_data=[]
    cloud_data=[]

    #-----------------------------------------------------------------------------------------------------------------
    # find the data and make it a 3D grid using mmr when weighting iportance
    # RH-data
    # use gas water burden to represent RH in 2D field
    
    if parameters[0]:
        # rh file name
        fname=os.path.join(path_to_folder,'aerocom3_'+model+'-met2010_AP3-CTRL_rh_ModelLevel_2010_monthly_T63L47.nc')
        # read the rh data (represented as gas water content burden) from the file to the variable
        rh_data = get_RH_burden(fname, fname_echam, fname_vphysc)
        
    else:
        # use SALSA base data
        fname=fname_ham
        # read the rh data (represented as gas water content burden) from the file to the variable
        rh_data = get_RH_burden(fname, fname_echam, fname_vphysc)
        
    #-----------------------------------------------------------------------------------------------------------------
    # mmr-data
    # use mass burden to represent mmr in 2D field
    
    if parameters[1]:
        # if mmr field is changed, use the model data to calculate mass burden
        mmr_file = os.path.join(path_to_folder,'aerocom3_'+model+'-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc')
        mmr_data = mass_burden(mmr_file, fname_vphysc)
    else:
        # otherwise use the salsa data to calculate mass burden
        mmr_file = os.path.join(path_to_folder,'aerocom3_ECHAM6.3-SALSA2.0-met2010_AP3-CTRL_mmrpm1_ModelLevel_2010_monthly_T63L47.nc')
        mmr_data = mass_burden(mmr_file, fname_vphysc)

    #-----------------------------------------------------------------------------------------------------------------
    # refractive indices data
    # files are saved in 4D to be available to use in other cases
    # all levels are identical as the indices and coefficients are same everywhere
    # therefore any level can be selected 
    
    if parameters[2]:
        # if refractive indices have been changed, use recalculated indices and coefficients
        # refractive indices consists of real part (refractive index) and imaginary part (absorption coefficient)
        
        # refractive indices file
        refrac_file = 'aerocom3_'+model+'-met2010_AP3-CTRL_refind_ModelLevel_2010_monthly_T63L47.nc'
        refrac_file = os.path.join(path_to_folder,refrac_file)
        
        # absorption coefficients file
        abs_file = 'aerocom3_'+model+'-met2010_AP3-CTRL_abs_ModelLevel_2010_monthly_T63L47.nc'
        abs_file = os.path.join(path_to_folder,abs_file)
        
        # If mmr is selected, then use pvols from that
        if parameters[1]:
            if not os.path.exists(refrac_file):
                create_refind(refrac_file,fname_tracer,fname_vphysc,True,model,mmr_file)   # create file containing refractive indices
            # get the 1st layer (all the layers are identical but some must be selected as the output must be 2D)
            refractive_ind_data = read_netcdf_variable(refrac_file,'ref_ind',0)

            if not os.path.exists(abs_file):
                create_refind(abs_file,fname_tracer,fname_vphysc,False,model,mmr_file,fname_ham)   # create file containing absorption coefficients
            # get the 1st layer from file (all the layers are identical but some must be selected as the output must be 2D)
            abs_data = read_netcdf_variable(abs_file,'abs_coef',0)
        
        # else just SALSA pvols
        else:
            # check if refractive indices and absorbtion coefficients are allready stored in a file
            if not os.path.exists(refrac_file):
                create_refind(refrac_file,fname_tracer,fname_vphysc,True,model)   # create file containing refractive indices
            # get the 1st layer (all the layers are identical but some must be selected as the output must be 2D)
            refractive_ind_data = read_netcdf_variable(refrac_file,'ref_ind',0)
        
            # same procedure for absorbtion coeficcients
            if not os.path.exists(abs_file):
                create_refind(abs_file,fname_tracer,fname_vphysc,False,model)   # create file containing absorption coefficients
            # get the 1st layer from file (all the layers are identical but some must be selected as the output must be 2D)
            abs_data = read_netcdf_variable(abs_file,'abs_coef',0)

    else:
        # use data from SALSA
        
        # naming the refractive indices file
        refrac_file = 'SALSA_refind_ModelLevel_2010_monthly_T63L47.nc'
        refrac_file = os.path.join(path_to_folder,refrac_file)
        
        # naming the absorption coefficients file
        abs_file = 'SALSA_abs_ModelLevel_2010_monthly_T63L47.nc'
        abs_file = os.path.join(path_to_folder,abs_file)
        
        # exactly the same steps for SALSA-based indices, except from the create_refind commands where the model name is missing and therefore
        if not os.path.exists(refrac_file):
            create_refind(refrac_file,fname_tracer,fname_vphysc,True)
        # read the data from file
        refractive_ind_data = read_netcdf_variable(refrac_file,'ref_ind',0)

        if not os.path.exists(abs_file):
            create_refind(abs_file,fname_tracer,fname_vphysc,False)
        abs_data = read_netcdf_variable(abs_file,'abs_coef',0)
    
    #-----------------------------------------------------------------------------------------------------------------
    # cloud fraction data
    # use maximum-random overlap method, the same as is used in DRE_offline_driver, to make into 2D
    
    if parameters[3]:
        # clt file name
        fname_clt=os.path.join(path_to_folder,'aerocom3_'+model+'-met2010_AP3-CTRL_clt_ModelLevel_2010_monthly_T63L47.nc')
        # read the clt data from the file to the variable
        with nc.Dataset(fname_clt, 'r', format='NETCDF4_CLASSIC') as cf:
            clt_var = cf.variables['clt']
            cloud_data = np.array(clt_var)
            # GISS and GFDL have cloud fractions as % so devide by 100 for that case
            if "GISS" in model or "GFDL" in model:
                cloud_data=cloud_data/100
            if cloud_data.ndim == 4:
                print("Detected 4D cloud fraction data. Estimating cover over all levels.")
                # add up cloud layers from bottom to top
                mr_clt = np.zeros_like(cloud_data[:, 0, :, :])  # (time, lat, lon)
                # use maximum-random overlap method for calculating total cover
                for l in reversed(range(47)):
                    cl = np.clip(cloud_data[:, l, :, :], 0, 1)
                    mr_clt = cl + (1 - cl) * mr_clt
                cloud_data = mr_clt
    else:
        # use SALSA base data
        with nc.Dataset(fname_echam, 'r', format='NETCDF4_CLASSIC') as sf:
            clt_var = sf.variables['aclcov']
            cloud_data = np.array(clt_var)


    #-----------------------------------------------------------------------------------------------------------------   

# when adding features, add here the additional feature data readings and remember to initialize a list for them and also add them to the dictionary
    
    # ensure data structure details
    print(f"RH data has shape of {rh_data.shape}")
    print(f"MMR data has shape of {mmr_data.shape}")
    print(f"Ref Ind data has shape of {refractive_ind_data.shape}")
    print(f"Abs Coef data has shape of {abs_data.shape}")
    print(f"Cloud data has shape of {cloud_data.shape}")

    # store the parameter data under corresponding variable
    root_data['rh']=rh_data
    root_data['mmr']=mmr_data
    root_data['refrac']=refractive_ind_data
    root_data['abs_coef']=abs_data
    root_data['clt']=cloud_data

    # check that values are not empty
    for key in root_data.keys():
        if root_data[key].size > 0:
            continue
        else:
            raise ValueError(f"Data under {key} is missing")
            
    return root_data
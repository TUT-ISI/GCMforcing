import numpy as np
import netCDF4 as nc
"""
by  Harri Kokkola / FMI

For forming and saving variabledata to netCDF4 files
Options for different data dimensions
"""

# Write a variable with rectilinear coordinates to NETCDF.
def write_2D_grid(file, var, lon, lat, varname=None):

    # open file for writing
    outdataset = nc.Dataset(file, 'w', format='NETCDF4_CLASSIC')
    
    # define longitudes and latitudes
    dlats  = outdataset.createDimension('lat', len(lat))
    dlons  = outdataset.createDimension('lon', len(lon))
    latitudes  = outdataset.createVariable('lat',np.float32, ('lat',))
    longitudes = outdataset.createVariable('lon',np.float32, ('lon',))

    # check if a dictionary variable is saved

    if isinstance(var,dict):

        ext=dict()
        for varname in var.keys():

            ext[varname] = outdataset.createVariable(varname,np.float64, ('lat','lon'))
            
            outdataset.variables['lat'].units = 'degrees_north'
            outdataset.variables['lon'].units = 'degrees_east'
            
            latitudes[:]=lat
            longitudes[:]=lon

        for varname in var.keys():

            ext[varname][:,:]=var[varname]

    # if only one variable is saved
    else:

        # if only one variable is saved

        varname=('varname' if varname is None else varname)

        ext    = outdataset.createVariable(varname,np.float64, ('lat','lon'))

        outdataset.variables['lat'].units = 'degrees_north'
        outdataset.variables['lon'].units = 'degrees_east'

        latitudes[:]=lat
        longitudes[:]=lon

        ext[:,:]=var

    outdataset.close()

    
def write_3D_grid(file, var, lon, lat, time, varname=None):

    # open file for writing
    outdataset = nc.Dataset(file, 'w', format='NETCDF4_CLASSIC')
    
    # define longitudes and latitudes
    dlats  = outdataset.createDimension('lat', len(lat))
    dlons  = outdataset.createDimension('lon', len(lon))
    dtimes = outdataset.createDimension('time', len(time))
    latitudes  = outdataset.createVariable('lat',np.float32, ('lat',))
    longitudes = outdataset.createVariable('lon',np.float32, ('lon',))
    times      = outdataset.createVariable('time', np.float32, ('time',))

    # check if a dictionary variable is saved

    if isinstance(var,dict):

        ext=dict()
        for varname in var.keys():

            ext[varname] = outdataset.createVariable(varname,np.float64, ('time','lat','lon'))
            
            outdataset.variables['lat'].units = 'degrees_north'
            outdataset.variables['lat'].standard_name = 'latitude'
            outdataset.variables['lat'].long_name = 'latitude'
            outdataset.variables['lat'].axis = 'Y'

            outdataset.variables['lon'].units = 'degrees_east'
            outdataset.variables['lon'].standard_name = 'longitude'
            outdataset.variables['lon'].long_name = 'longitude'
            outdataset.variables['lon'].axis = 'X'
           
            outdataset.variables['time'].units = '1'
            outdataset.variables['time'].units = 'days since 2009-01-01 00:00:00'
            outdataset.variables['time'].standard_name = 'time'
            outdataset.variables['time'].calendar = 'proleptic_gregorian'
            
            latitudes[:]=lat
            longitudes[:]=lon
            times[:]=time

        for varname in var.keys():

            ext[varname][:,:,:]=var[varname]

    # if only one variable is saved
    else:

        # if only one variable is saved

        varname=('varname' if varname is None else varname)

        ext    = outdataset.createVariable(varname,np.float64, ('time','lat','lon'))

        outdataset.variables['lat'].units = 'degrees_north'
        outdataset.variables['lat'].standard_name = 'latitude'
        outdataset.variables['lat'].long_name = 'latitude'
        outdataset.variables['lat'].axis = 'Y'
 
        outdataset.variables['lon'].units = 'degrees_east'
        outdataset.variables['lon'].standard_name = 'longitude'
        outdataset.variables['lon'].long_name = 'longitude'
        outdataset.variables['lon'].axis = 'X'

        outdataset.variables['time'].units = 'days since 2009-01-01 00:00:00'
        outdataset.variables['time'].standard_name = 'time'
        outdataset.variables['time'].calendar = 'proleptic_gregorian'
        
        latitudes[:]=lat
        longitudes[:]=lon
        times[:]=time

        ext[:,:,:]=var

    outdataset.close()

def write_4D_grid(file, var, lon, lat, lev, time, varname=None):

    # open file for writing
    outdataset = nc.Dataset(file, 'w', format='NETCDF4_CLASSIC')
    
    # define longitudes and latitudes
    dlats  = outdataset.createDimension('lat', len(lat))
    dlons  = outdataset.createDimension('lon', len(lon))
    dlevs  = outdataset.createDimension('lev', len(lev))
    dtimes = outdataset.createDimension('time', len(time))
    latitudes  = outdataset.createVariable('lat',np.float32, ('lat',))
    longitudes = outdataset.createVariable('lon',np.float32, ('lon',))
    levels     = outdataset.createVariable('lev',np.float32, ('lev',))
    times      = outdataset.createVariable('time', np.float32, ('time',))

    # check if a dictionary variable is saved

    if isinstance(var,dict):

        ext=dict()
        for varname in var.keys():

            ext[varname] = outdataset.createVariable(varname,np.float64, ('time','lev','lat','lon'))
            
            outdataset.variables['lat'].units = 'degrees_north'
            outdataset.variables['lat'].standard_name = 'latitude'
            outdataset.variables['lat'].long_name = 'latitude'
            outdataset.variables['lat'].axis = 'Y'

            outdataset.variables['lon'].units = 'degrees_east'
            outdataset.variables['lon'].standard_name = 'longitude'
            outdataset.variables['lon'].long_name = 'longitude'
            outdataset.variables['lon'].axis = 'X'

            outdataset.variables['lev'].units = 'level'
            outdataset.variables['lev'].standard_name = 'hybrid_sigma_pressure'
            outdataset.variables['lev'].long_name = 'hybrid level at layer midpoints'
            outdataset.variables['lev'].units = 'level'
            outdataset.variables['lev'].positive = 'down'
            outdataset.variables['lev'].formula = 'hyam hybm (mlev=hyam+hybm*aps)'
            outdataset.variables['lev'].formula_terms = 'ap: hyam b: hybm ps: aps'
            
            outdataset.variables['time'].units = '1'
            outdataset.variables['time'].units = 'days since 2009-01-01 00:00:00'
            outdataset.variables['time'].standard_name = 'time'
            outdataset.variables['time'].calendar = 'proleptic_gregorian'
            
            latitudes[:]=lat
            longitudes[:]=lon
            levels[:]=lev
            times[:]=time

        for varname in var.keys():

            ext[varname][:,:,:,:]=var[varname]

    # if only one variable is saved
    else:

        # if only one variable is saved

        varname=('varname' if varname is None else varname)

        ext    = outdataset.createVariable(varname,np.float64, ('time','lev','lat','lon'))

        outdataset.variables['lat'].units = 'degrees_north'
        outdataset.variables['lat'].standard_name = 'latitude'
        outdataset.variables['lat'].long_name = 'latitude'
        outdataset.variables['lat'].axis = 'Y'
 
        outdataset.variables['lon'].units = 'degrees_east'
        outdataset.variables['lon'].standard_name = 'longitude'
        outdataset.variables['lon'].long_name = 'longitude'
        outdataset.variables['lon'].axis = 'X'

        #outdataset.variables['lev'].units = 'level'
        outdataset.variables['lev'].units = '1'
        outdataset.variables['lev'].standard_name = 'hybrid_sigma_pressure'
        outdataset.variables['lev'].long_name = 'hybrid level at layer midpoints'
        outdataset.variables['lev'].units = 'level'
        outdataset.variables['lev'].positive = 'down'
        outdataset.variables['lev'].formula = 'hyam hybm (mlev=hyam+hybm*aps)'
        outdataset.variables['lev'].formula_terms = 'ap: hyam b: hybm ps: aps'

        outdataset.variables['time'].units = 'days since 2009-01-01 00:00:00'
        outdataset.variables['time'].standard_name = 'time'
        outdataset.variables['time'].calendar = 'proleptic_gregorian'
        
        latitudes[:]=lat
        longitudes[:]=lon
        levels[:]=lev
        times[:]=time

        ext[:,:,:,:]=var

    outdataset.close()

    
    

import csv
import numpy as np
import netCDF4 as nc
import calendar
import sys, os
from calendar import monthrange
from netCDF4 import date2num
from datetime import datetime, timedelta, date
from cdo import Cdo as CDO
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

"""
by  Harri Kokkola / FMI 
    Atte Laakso / Aalto University

Regrids given (3D) file based on salsagrid.nc
"""

os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
cdo_path = os.getenv('CDO')
cdo = CDO()

# Write a variable with rectilinear coordinates to NETCDF (3D: time, lat, lon)
def write_3D_grid(file, var, lon, lat, time, varname, unit, long_name):

    # open file for writing
    outdataset = nc.Dataset(file, 'w', format='NETCDF4_CLASSIC')
    
    # define dimensions
    dlats  = outdataset.createDimension('lat', len(lat))
    dlons  = outdataset.createDimension('lon', len(lon))
    dtimes = outdataset.createDimension('time', len(time))
    
    # create variables
    latitudes  = outdataset.createVariable('lat', np.float32, ('lat',))
    longitudes = outdataset.createVariable('lon', np.float32, ('lon',))
    times      = outdataset.createVariable('time', np.float32, ('time',))
    
    varname = ('varname' if varname is None else varname)
    ext = outdataset.createVariable(varname, np.float64, ('time', 'lat', 'lon'), zlib=True)
    
    # variable attributes
    outdataset.variables['lat'].units = 'degrees_north'
    outdataset.variables['lat'].standard_name = 'latitude'
    outdataset.variables['lat'].long_name = 'latitude'
    outdataset.variables['lat'].axis = 'Y'
    
    outdataset.variables['lon'].units = 'degrees_east'
    outdataset.variables['lon'].standard_name = 'longitude'
    outdataset.variables['lon'].long_name = 'longitude'
    outdataset.variables['lon'].axis = 'X'
    
    outdataset.variables['time'].units = 'days since 2001-01-01'
    outdataset.variables['time'].standard_name = 'time'
    outdataset.variables['time'].calendar = 'proleptic_gregorian'

    outdataset.variables[varname].units = unit
    outdataset.variables[varname].long_name = long_name

    # assign data
    latitudes[:] = lat
    longitudes[:] = lon
    times[:] = time
    print(var.shape, time.shape, lat.shape, lon.shape)
    ext[:, :, :] = var
    
    outdataset.close()

# Model name
model = sys.argv[1]

# Variable to be processed
variable = sys.argv[2]

# Folder that contains this file
p = sys.argv[3]

# Construct input filename
input_fname = f'aerocom3_{model}-met2010_AP3-CTRL_{variable}_ModelLevel_2010_monthly.nc'
input_fname = os.path.join(p, input_fname)

output_fname = f'aerocom3_{model}-met2010_AP3-CTRL_{variable}_ModelLevel_2010_monthly_T63L47.nc'
output_fname = os.path.join(p, output_fname)

modeldata = nc.Dataset(input_fname, 'r', format='NETCDF4_CLASSIC')
print(modeldata.variables[variable].units)

# Get long_name and unit of the variable
if hasattr(modeldata.variables[variable], 'long_name'):
    long_name = modeldata.variables[variable].long_name
    unit = modeldata.variables[variable].units
else:
    long_name = variable
    unit = ''

# Read grid parameters
lon = modeldata.variables['lon'][:]
lat = modeldata.variables['lat'][:]
time = modeldata.variables['time'][:]

variable_value = modeldata.variables[variable][:]  # expected shape: (time, lat, lon)

# Interpolation on horizontal grid if needed

# Load SALSA grid (assuming salsa grid is 2D lat, lon)
salsadata = nc.Dataset('salsagrid.nc', 'r', format='NETCDF_CLASSIC')
lat_salsa = salsadata.variables['lat'][:]
lon_salsa = salsadata.variables['lon'][:]

# Prepare output array
z_new = np.zeros((len(time), len(lat_salsa), len(lon_salsa)))

# Flatten input grid for interpolation
points = np.array([(lo, la) for la in lat for lo in lon])  # mesh grid flattening
for t in range(len(time)):
    values = variable_value[t, :, :].flatten()
    
    # Create target grid points
    lon_salsa_2d, lat_salsa_2d = np.meshgrid(lon_salsa, lat_salsa)
    target_points = np.array([lon_salsa_2d.flatten(), lat_salsa_2d.flatten()]).T
    
    # Interpolate using griddata (cubic) and interpolate any nan values using (nearest)
    interp_values = griddata(points, values, target_points, method='cubic')
    nan_mask = np.isnan(interp_values)
    if np.any(nan_mask):
        interp_values[nan_mask] = griddata(points, values, target_points[nan_mask], method='nearest')
    
    # Reshape to 2D
    z_new[t, :, :] = interp_values.reshape(len(lat_salsa), len(lon_salsa))

print("New shape", z_new.shape)

write_3D_grid('tmp2.nc', z_new, lon_salsa, lat_salsa, time, variable, unit, long_name)

# Use CDO to remap conservatively (or just copy if remap is not needed)
cdo.copy(input='-remapcon,salsagrid.nc tmp2.nc', output=output_fname)

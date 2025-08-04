import csv
import numpy as np
import netCDF4 as nc
import calendar
import sys, os
from calendar import monthrange
from netCDF4 import date2num
from datetime import datetime, timedelta, date
from cdo import Cdo as CDO
from scipy import interpolate
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

"""
by  Harri Kokkola / FMI 

Regrids given file based on salsagrid.nc
"""

# CDO initialization
os.environ['CDO']='/appl/spack/v018/install-tree/gcc-8.5.0/cdo-2.0.5-zpo6xz/bin/cdo' 
cdo_path=os.getenv('CDO')
cdo = CDO()

# Write a variable with rectilinear coordinates to NETCDF.
def write_4D_grid(file, var, lon, lat, lev, time, hyai, hybi, hyam, hybm, varname, unit, long_name):

    
    # open file for writing
    outdataset = nc.Dataset(file, 'w', format='NETCDF4_CLASSIC')
    
    # define longitudes and latitudes
    dlats  = outdataset.createDimension('lat',len(lat))
    dlons  = outdataset.createDimension('lon',len(lon))
    dlevs  = outdataset.createDimension('lev',len(lev))
    dtimes = outdataset.createDimension('time',len(time))
    dhyai  = outdataset.createDimension('hyai',len(hyai))
    dhybi  = outdataset.createDimension('hybi',len(hybi))
    dhyam  = outdataset.createDimension('hyam',len(hyam))
    dhybm  = outdataset.createDimension('hybm',len(hybm))
    latitudes  = outdataset.createVariable('lat',np.float32, ('lat',))
    longitudes = outdataset.createVariable('lon',np.float32, ('lon',))
    levels     = outdataset.createVariable('lev',np.float32, ('lev',))
    times      = outdataset.createVariable('time', np.float32, ('time',))
    hybridais =  outdataset.createVariable('hyai', np.float32, ('hyai',))
    hybridbis =  outdataset.createVariable('hybi', np.float32, ('hybi',))
    hybridams =  outdataset.createVariable('hyam', np.float32, ('hyam',))
    hybridbms =  outdataset.createVariable('hybm', np.float32, ('hybm',))
    # check if a dictionary variable is saved

    varname=('varname' if varname is None else varname)
    
    ext    = outdataset.createVariable(varname,np.float64, ('time','lev','lat','lon'),zlib=True)
    
    outdataset.variables['lat'].units = 'degrees_north'
    outdataset.variables['lat'].standard_name = 'latitude'
    outdataset.variables['lat'].long_name = 'latitude'
    outdataset.variables['lat'].axis = 'Y'
    
    outdataset.variables['lon'].units = 'degrees_east'
    outdataset.variables['lon'].standard_name = 'longitude'
    outdataset.variables['lon'].long_name = 'longitude'
    outdataset.variables['lon'].axis = 'X'
    
    outdataset.variables['lev'].standard_name = "hybrid_sigma_pressure" ;
    outdataset.variables['lev'].long_name = "hybrid level at layer midpoints" ;
    outdataset.variables['lev'].formula = "hyam hybm (mlev=hyam+hybm*aps)" ;
    outdataset.variables['lev'].formula_terms = "ap: hyam b: hybm ps: aps" ;
    outdataset.variables['lev'].units = "level" ;
    outdataset.variables['lev'].positive = "down" ;

    outdataset.variables['hyai'].long_name = "hybrid A coefficient at layer interfaces" 
    outdataset.variables['hyai'].units = "Pa" 

    outdataset.variables['hybi'].long_name = "hybrid B coefficient at layer interfaces" 
    outdataset.variables['hybi'].units = "1" 
    
    outdataset.variables['hyam'].long_name = "hybrid A coefficient at layer midpoints" 
    outdataset.variables['hyam'].units = "Pa" 

    outdataset.variables['hybm'].long_name = "hybrid B coefficient at layer midpoints" 
    outdataset.variables['hybm'].units = "1" 
   
    outdataset.variables['time'].units = 'day as %Y-%m-%d'
    outdataset.variables['time'].units = 'days since 2001-01-01'
    outdataset.variables['time'].standard_name = 'time'
    outdataset.variables['time'].calendar = 'proleptic_gregorian'

    outdataset.variables[varname].units = unit
    outdataset.variables[varname].long_name = long_name

    latitudes[:]=lat
    longitudes[:]=lon
    levels[:]=lev
    times[:]=time
    hybridais[:] = hyai
    hybridbis[:] = hybi
    hybridams[:] = hyam
    hybridbms[:] = hybm
    print( var.shape, time.shape, lev.shape, lat.shape, lon.shape)
    ext[:,:,:,:]=var
    

    outdataset.close()

# Model name
model=sys.argv[1]

# Variable to be processed
variable=sys.argv[2]

# Folder that contains this file
p = sys.argv[3]

# Construct input filename
input_fname='aerocom3_'+model+'-met2010_AP3-CTRL_'+variable+'_ModelLevel_2010_monthly.nc'
input_fname=os.path.join(p,input_fname)

output_fname='aerocom3_'+model+'-met2010_AP3-CTRL_'+variable+'_ModelLevel_2010_monthly_T63L47.nc'
output_fname=os.path.join(p,output_fname)

modeldata=nc.Dataset(input_fname, 'r', format='NETCDF4_CLASSIC')    
print(modeldata.variables[variable].units)

# Get long_name and the unit of the variable
# First check if it exists
if hasattr(modeldata.variables[variable], 'long_name'):
    long_name=modeldata.variables[variable].long_name
    unit=modeldata.variables[variable].units
else:
# If not, use empty string for the unit and the variable name for long name
    long_name=variable
    unit=''

lev_str='lev'

if model=='INCA':
    lev_str='pres'

# Read in the model grid parameters
lon=modeldata.variables['lon'][:]
lat=modeldata.variables['lat'][:]
lev_tmp=modeldata.variables[lev_str][:]

lev=range(1,len(lev_tmp)+1)
variable_value2=modeldata.variables[variable][:]

# reverse levels for certain models
if model=='OsloCTM3v1.01' or model=='GFDL-AM4' or model=='GISS-ModelE2p1p1-OMA' or model=='INCA' or model=='MIROC-SPRINTARS':
    z=variable_value2[:,::-1,:,:]
    print( 'reverse', model)
if model=='CAM5-ATRAS':
    z=variable_value2
    print( 'forward', model)


# Get the SALSA grid parameters
salsadata=nc.Dataset('salsagrid.nc', 'r', format='NETCDF_CLASSIC')
lev_salsa=salsadata.variables['lev'][:]
hyai=salsadata.variables['hyai'][:]
hybi=salsadata.variables['hybi'][:]
hyam=salsadata.variables['hyam'][:]
hybm=salsadata.variables['hybm'][:]
time=modeldata.variables['time'][:]

L47=np.linspace(1.0,len(lev_tmp),num=47)
L60=np.linspace(1.0,len(lev_tmp),num=60)
z_new=np.zeros((len(time),len(L47),len(lat),len(lon)))
for i in range(0,len(lon)):

    for j in range(0,len(lat)):

        for k in range(0,len(time)):

            variable_value=np.squeeze(z[k,:,j,i])

            x=np.squeeze(L47)
            xi=np.squeeze(lev)
            y=np.squeeze(lat)
            yi=y

            z_new[k,:,j,i]=griddata(xi, variable_value.flatten(), x, method='cubic')
print("New shape", z_new.shape)

write_4D_grid('tmp2.nc', z_new, lon, lat, lev_salsa, time, hyai, hybi, hyam, hybm, variable, unit, long_name)

cdo.copy(input='-remapcon,salsagrid.nc tmp2.nc', output=output_fname)

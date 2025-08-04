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
cdo = CDO()

# Write a variable with rectilinear coordinates to NETCDF.
def write_4D_grid(file, var, aps, lon, lat, lev, time, hyai, hybi, hyam, hybm, varname, unit, long_name):

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
    psurf  = outdataset.createVariable('aps',np.float64, ('time','lat','lon'))
    
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

    outdataset.variables['aps'].units = 'Pa'
    outdataset.variables['aps'].standard_name = 'surface_air_pressure'
    
    latitudes[:]=lat
    longitudes[:]=lon
    levels[:]=lev
    times[:]=time
    print len(hyai), len(hyam)
    hybridais[:] = hyai
    hybridbis[:] = hybi
    hybridams[:] = hyam
    hybridbms[:] = hybm
    print aps.shape
    psurf[:,:,:] = aps
    ext[:,:,:,:]=var
    

    outdataset.close()


# tracer file
input_fname=sys.argv[1]

# vphysc file
output_fname=sys.argv[2]

variable=sys.argv[3]

unit=sys.argv[4]

long_name=sys.argv[5]

modeldata=nc.Dataset(input_fname, 'r', format='NETCDF4_CLASSIC')    
salsadata=nc.Dataset('salsagrid.nc', 'r', format='NETCDF_CLASSIC')
lon=modeldata.variables['lon'][:]
lat=modeldata.variables['lat'][:]
lev_salsa=salsadata.variables['lev'][:]
lev_tmp=modeldata.variables['lev'][:]
lev=range(1,len(lev_tmp)+1)
hyai=modeldata.variables['a_bnds'][:]
hybi=modeldata.variables['b_bnds'][:]
hyam2=modeldata.variables['a'][:]
hybm2=modeldata.variables['b'][:]
aps2=modeldata.variables['ps'][:]
hyai2=[]
hybi2=[]
for i in range(0,len(hyai)):
    hyai2.append(hyai[i][0])
    hybi2.append(hybi[i][0])
hyai2.append(hyai[i,1])
hybi2.append(hybi[i,1])

hyai=hyai2[::-1]
hybi=hyai2[::-1]
hyam=hyam2[::-1]
hybm=hybm2[::-1]
aps=aps2[::-1]

time=modeldata.variables['time'][:]

variable_value2=modeldata.variables[variable][:]

variable_value=variable_value2[:,::-1,:,:]

L47=np.linspace(1.0,60.0,num=47)

print variable_value.shape, time.shape
#print lev.shape, lat.shape, lon.shape

variable_final=griddata((time.flatten(), lev_salsa, lat, lon), variable_value, (time, lev, lat, lon), method='cubic')


#write_4D_grid('tmp2.nc', variable_value, aps, lon, lat, lev, time, hyai, hybi, hyam, hybm, variable, unit, long_name)
#cdo.copy(input='-setzaxis,vct_L60.txt test.nc', output='tmp2.nc')
#cdo.copy(input='-remapeta,vct_L47.txt tmp2.nc', output=output_fname)


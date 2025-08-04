#!/usb/bin/python3
# Written for Python v2.7.6
#    Harri Kokkola / FMI 
"""

Calculates the relative wet radius above both land and ocean.
 
"""
from netCDF4 import Dataset
from dpwet import wet_diameter
from aod import calculate_aod
#from map_visuals import plot_map
import sys, os
import numpy as np
import netCDF4 as nc
from read_tracers import read_aero_bins
from refractive_index import *
from write_netcdf import *
from cdo import Cdo as CDO
from cdnc import cloud_activation
cdo = CDO()
print(len(sys.argv))
# Check if the command has the correct amount of input files
if len(sys.argv) != 11:
    print("\nInvalid number of inputs.\n   Usage: extinction_driver.py tracer_file.nc vphysc_file.nc ham_file.nc activ_file.nc echam_file.nc extinction_output.nc absorption_output.nc activation_output.nc wavelength relative_humidity \n\n")
    exit()

# Filenames read from command line input:

# tracer file
fname_tracer=sys.argv[1]

# vphysc file
fname_vphysc=sys.argv[2]

# ham_file
fname_ham=sys.argv[3]

# activ_file
fname_activ=sys.argv[4]

# echam_file
fname_echam=sys.argv[5]

cdo.copy(
    input = ' '.join([
        '-sp2gp',
        '-selname,st',
        fname_echam,
    ]),
    output = 'tmp_temp.nc',
)

# extinction output file
fname3=sys.argv[6]
# absorption output file
fname4=sys.argv[7]
# cloud activation output file
fname5=sys.argv[8]

# wavelength in nm
lmbd_str=sys.argv[9]
wavelength=int(lmbd_str)*1e-9

if(sys.argv[10] =='ambient'):

    relhum_fixed=-999

else:
    
    relhum_fixed=int(sys.argv[10])/100.0

model='SALSA'

def main():

    # string for a fixed rh (e.g. rh80 for 80%)
    rh_str='rh'+str(int(relhum_fixed*100))

    # if RH=0, hydration is not calculated
    if relhum_fixed==0:
        rh_str='dry'

    # if ambient model RH is used
    if relhum_fixed<0:
        rh_str=''

    # Check if files exist 
    if not os.path.exists(fname_tracer):
        print("Tracer file not found.")        
        exit()

    if not os.path.exists(fname_vphysc):
        print("vphysc file not found.")
        exit()

    # Read in relative humidity
        
    griddata=nc.Dataset(fname_vphysc, 'r', format='NETCDF4_CLASSIC')    
    lon=griddata.variables['lon'][:]
    lat=griddata.variables['lat'][:]
    lev=griddata.variables['lev'][:]
    time=griddata.variables['time'][:]

    griddim=3

    tracerfile=fname_tracer
    vphyscfile=fname_vphysc

    # get number and volume concentrations of individual species
    print('Reading in number and volume concentrations')
    znaero,zvols,time = read_aero_bins(tracerfile,vphyscfile,lon,lat)

    # get grid cell height for calculating AOD
    print('Reading in grid box height from vphysc_file.nc')
    vphyscdata = Dataset(vphyscfile, 'r')
    h = vphyscdata.variables['grheightm1']
    zgrheight = h

    # calculate full-level pressure (pfull) from half-level values (phalf)
    phalf = vphyscdata.variables['aphm1']
    pfull = (phalf[:,0:47,:,:] + phalf[:,1:48,:,:])/2.
    zapm1 = pfull
    
    # get relative humidity for model levels
    print('Reading in relative humidity from ham_file.nc')
    hamdata=nc.Dataset(fname_ham, 'r', format='NETCDF4_CLASSIC')    
    relhum = hamdata.variables['relhum']
    zrh = relhum

    # get updraft velocities for model levels
    print('Reading updraft velocities from activ_file.nc')
    activdata=nc.Dataset(fname_activ, 'r', format='NETCDF4_CLASSIC')
    vervel = activdata.variables['W']
    zw = vervel

    # get specific humidity for model levels
    echamdata=nc.Dataset(fname_echam, 'r', format='NETCDF_CLASSIC')
    specific_humidity=echamdata.variables['q']
    zqm1 = specific_humidity
    
    
    # get temperature for model levels
    temperaturedata=nc.Dataset('tmp_temp.nc', 'r', format='NETCDF_CLASSIC')
    temperature=temperaturedata.variables['st']
    ztm1 = np.array(temperature)

    os.remove('tmp_temp.nc')
    
    # calculate the wet diameter of particles in each size bin
    print('Calculating the wet diameter')
    dwet,ddry=wet_diameter(znaero,zvols,zrh,model,relhum_fixed)

    # calculate the refractive index based on volume concentrations
    print('Calculating refractive indices')
    #nr,ni=refractive_index(zvols)

    print('Calculating AOD')
    lut_file='lut_optical_properties.nc'
    #aerosol_optical_depth,extinction,absorption=calculate_aod(nr,ni,dwet,znaero,zgrheight,lut_file, model, wavelength)

    print('Calculating number of activated droplets')
    cdnc=cloud_activation(znaero, zvols, ztm1, zapm1, zqm1, zw, 1)#, 1)
    
    print('Saving data to NetCDF files')

    if(griddim == 2):
        
        write_2D_grid('dwet.nc',dwet,lon,lat)

        write_2D_grid('ddry.nc',ddry,lon,lat)
        
        write_2D_grid('nr.nc',nr,lon,lat)
        
        write_2D_grid('ni.nc',ni,lon,lat)
        
        write_2D_grid('rh.nc',rh,lon,lat)
        
        write_2D_grid(fname3,aerosol_optical_depth,lon,lat,'AOD')

    elif(griddim == 3):

	#write_4D_grid('dwet.nc',dwet,lon,lat,lev,time,'dpwet')
        #write_4D_grid(fname3,extinction,lon,lat,lev,time,'ec'+lmbd_str+rh_str+'aer')
        #write_4D_grid(fname4,absorption,lon,lat,lev,time,'ac'+lmbd_str+rh_str+'aer')
        write_4D_grid(fname5,cdnc,lon,lat,lev,time,'CDNC')

# do not run main if imported
if __name__ == "__main__":
    main()
        

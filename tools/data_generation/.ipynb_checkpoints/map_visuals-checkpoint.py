# Written for Python v2.7.6
# by Tomi Laaksoviita / FMI
#    Harri Kokkola / FMI 
#import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pylab import *

def plot_map(zparam,lon,lat,zlim=None,ncolors=None):

    # If z-scale is not defined, use min and max values of the parameter for the color scale
    zlim=([np.amin(zparam),np.amax(zparam)] if zlim is None else zlim)

    print(zlim)
    # Figure size
    fig1 = plt.figure(1,figsize=(12,4))
    ax1 = fig1.add_axes([0.1,0.1,0.6,0.8])

    bmap = bm.Basemap(llcrnrlon=-180.,llcrnrlat=-80,urcrnrlon=180.,urcrnrlat=80,projection='mill',ax=ax1,lon_0=0)
    bmap.drawcoastlines(linewidth=1.)
    bmap.drawparallels(np.arange(-80.,80.,20.),labels=[1,0,0,0])
    bmap.drawmeridians(np.arange(-180.,240.,60.),labels=[0,0,0,1])
    zparamsh, lonsh = bm.shiftgrid(180.,zparam,lon)

    # shift the map for a more "European" look
    lonsh = lonsh - 360.
    glons,glats = np.meshgrid(lonsh,lat)
    xx,yy = bmap(glons,glats)
    
    # calculate ticks based on number of colors required
    ncolors=(5 if ncolors is None else ncolors)
    zstep=(zlim[1]-zlim[0])/ncolors
    
    # plot map
    pl1 = bmap.contourf(xx,yy,zparamsh,np.arange(zlim[0],zlim[1],zstep),cmap=plt.get_cmap('YlOrRd'),extend='both')
             
    # Make colorbar
    cb=bmap.colorbar(pl1, "right", size="5%", pad="10%")
    cb.set_ticks(np.arange(zlim[0],zlim[1],zstep))
    
    plt.show()


 

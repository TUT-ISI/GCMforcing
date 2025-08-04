import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import matplotlib.pyplot as plt
import sys, os
import netCDF4 as nc
import matplotlib.colors as mcolors

"""
Written for Python 3
by Tomi Laaksoviita / FMI
   Harri Kokkola / FMI
 
edited by Atte Laakso / Aalto University 
    to use Cartopy instead of Basemap and to give colormap as argument
    to use discrete data when label states to have FI-data in zparam
    to use clipped continous map for R2 and RMSE plots
"""

def plot_map(zparam,lon,lat,outfile,label,colormap,feature_names,zlim=None,ncolors=None):
    is_feature_plot = 'feature' in label.lower()

    # Prepare feature importance map or metric map
    if is_feature_plot:
        zparam = np.round(zparam).astype(int)  # Convert to discrete integers
        unique_vals = np.unique(zparam)
    else:
        zparam = np.clip(zparam, 0, None)
        unique_vals = None

    # Setup limits for maps
    zlim = [np.amin(zparam), np.amax(zparam)] if zlim is None else zlim
    ncolors = 4 if ncolors is None else ncolors

    # Plot general setup
    fig, ax = plt.subplots(figsize=(12, 5), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_global()
    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False

    # Meshgrid
    glons, glats = np.meshgrid(lon, lat)

    # Plotting values
    if is_feature_plot:
        # Get unique values and total number of features in data
        unique_vals = np.unique(zparam)
        n_features = len(unique_vals)

        # Create boundaries
        boundaries = np.arange(np.min(unique_vals) - 0.5, np.max(unique_vals) + 1.5, 1)
        n_bins = len(boundaries) - 1

        # Ensure enough colors
        ncolors = max(ncolors, n_bins)

        # Create a discrete colormap with enough colors
        cmap = plt.get_cmap(colormap, ncolors)
        norm = mcolors.BoundaryNorm(boundaries, cmap.N)

        # Plot with discrete colors
        cf = ax.pcolormesh(glons, glats, zparam, cmap=cmap, norm=norm, shading='auto',
                        transform=ccrs.PlateCarree())

        # Set ticks at actual feature indices
        cb_ticks = unique_vals
        cb_labels = [str(i) if feature_names is None else feature_names[i] for i in unique_vals]

    else:
        vmax = (15 if 'RMSE' in label else 1) # Setting consistent limits for RMSE and R²
        cf = ax.pcolormesh(glons, glats, zparam, cmap=colormap, shading='auto',
                           transform=ccrs.PlateCarree(), vmin=0, vmax=vmax)

        cb_ticks = np.linspace(0, vmax, 6)
        cb_labels = None

    # Add colorbar
    cbar = plt.colorbar(cf, ax=ax, orientation='vertical', ticks=cb_ticks, shrink=0.6, pad=0.05)
    cbar.set_label(label)
    if cb_labels:
        cbar.set_ticks(cb_ticks)
        cbar.set_ticklabels(cb_labels)

    plt.savefig(outfile, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved to: {outfile}")
    return
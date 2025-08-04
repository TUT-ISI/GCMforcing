# +
import numpy as np

"""
by  Atte Laakso / Aalto University

Functions for combining 4D field's all vertical levels into one
"""

def combine_aods(aods,lat,lon,time):
    print("Combining AOD levels")
    # Initialize the combined array with shape (time, lat, lon)
    combined_aods = np.zeros((len(time), len(lat), len(lon)))
    # combine levels for every point in time and space
    for t in range(len(time)):
        for i in range(len(lat)):
            for j in range(len(lon)):
                    combined_aods[t, i, j] = np.sum(aods[t, :, i, j])
                    
    return combined_aods

def combine_ext(extinction,lat,lon,time):
    print("Combining extinction levels")
    # Initialize the combined array with shape (time, lat, lon)
    combined_ext = np.zeros((len(time), len(lat), len(lon)))
    # combine levels for every point in time and space
    for t in range(len(time)):
        for i in range(len(lat)):
            for j in range(len(lon)):
                    combined_ext[t, i, j] = np.sum(extinction[t, :, i, j])
                    
    return combined_ext

def combine_abs(absorption,lat,lon,time):
    print("Combining absorption levels")
    # Initialize the combined array with shape (time, lat, lon)
    combined_abs = np.zeros((len(time), len(lat), len(lon)))
    # combine levels for every point in time and space
    for t in range(len(time)):
        for i in range(len(lat)):
            for j in range(len(lon)):
                    combined_abs[t, i, j] = np.sum(absorption[t, :, i, j])
                    
    return combined_abs
# -



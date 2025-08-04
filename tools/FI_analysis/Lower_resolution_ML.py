import os, sys
import netCDF4 as nc
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import time
from ML_additionals import *
import dask
from dask.distributed import wait
from write_netcdf import *
from visualizations import plot_map

def low_res_ML(x_train, y_train, x_test, y_test, original_shape, domain_lat, domain_lon,
                fnames, lon, lat, use_log, rs, modelsel,
               output_dir='FI_results', client=None):
    """
    by  Atte Laakso / Aalto University

    Train given model on domain-wise subsets of gridded data,
    calculate feature importances, R2 scores and rmse, then save these results.
    
    Parameters:
    - x_train, y_train, x_test, y_test: pandas DataFrames or arrays
    - original_shape: tuple of (lat_size, lon_size) original grid dimensions
    - domain_lat, domain_lon: integers, size of subdomains in lat and lon
    - client: Dask client object
    - fnames: list of feature names
    - lon, lat: coordinate arrays
    - use_log: boolean, whether data is log transformed
    - output_dir: string, directory for saving results
    - modelsel: string, defines model to be used (see code below + ML_basis)
    
    Returns: None but saves the calculated data to outputfiles
    """

    # keep track of time
    start_time = time.time()

    print(f"Defining {modelsel} models for domains")
    # form gridded data for every set
    x_test_domains = form_gridwise_data( x_test , original_shape , domain_lat , domain_lon )
    y_test_domains = form_gridwise_data( y_test , original_shape , domain_lat , domain_lon )
    x_train_domains = form_gridwise_data( x_train , original_shape , domain_lat , domain_lon )
    y_train_domains = form_gridwise_data( y_train , original_shape , domain_lat , domain_lon )
    
    print("Calculating feature importances")
    # define number of domains
    num_domains_x = original_shape[0]//domain_lat
    num_domains_y = original_shape[1]//domain_lon
    # define ranges over the original dimensions to form the domainwise ML models
    I = range(num_domains_x)
    J = range(num_domains_y)
    
    # distribute tasks
    tasks           = []   # futures
    domain_indices  = []   # (i, j)-key for each future

    for i in range(num_domains_x):
        for j in range(num_domains_y):
            xt, yt, xte, yte = x_train_domains[i][j], y_train_domains[i][j], \
                            x_test_domains[i][j], y_test_domains[i][j]

            future = client.submit(get_FI, xt, yt, xte, yte,
                                fnames, rs, modelsel)
            tasks.append(future)
            domain_indices.append((i, j))

    # gather results
    results = client.gather(tasks)

    # free up some memory
    del x_train, x_train_domains, y_train, y_train_domains, x_test, x_test_domains, y_test, y_test_domains

    # map fi, r2 and rmse:s back to the grid
    fi_scores, r2_scores, rmse_scores = {}, {}, {}
    for (i, j), (fi, r2, rmse_val) in zip(domain_indices, results):
        fi_scores[(i, j)]   = fi
        r2_scores[(i, j)]   = r2
        rmse_scores[(i, j)] = rmse_val
        
    # print out means
    print(f"Mean R2 score is {np.mean(list(r2_scores.values())):.4f}")
    print(f"Mean RMSE is {np.mean(list(rmse_scores.values())):.4f}")

    # dictionaries to arrays
    r2_array = np.array([[r2_scores[(i, j)] for j in range(num_domains_y)] for i in range(num_domains_x)])
    rmse_array = np.array([[rmse_scores[(i, j)] for j in range(num_domains_y)] for i in range(num_domains_x)])

    # FI scores array
    fi_scores_array = np.empty((num_domains_x, num_domains_y), dtype=object)
    for i in range(num_domains_x):
        for j in range(num_domains_y):
            fi_scores_array[i, j] = fi_scores[(i, j)]

    # initialize FI index array
    FI_indices = np.full((num_domains_x, num_domains_y, len(fnames)), -1, dtype=int)

    print("Processing feature importances per domain")

    for i in range(num_domains_x):
        for j in range(num_domains_y):
            keys = list(fi_scores_array[i, j].keys())
            values = list(fi_scores_array[i, j].values())

            # sort features by importance descending
            sorted_indices = sorted(range(len(values)), key=lambda idx: values[idx], reverse=True)

            for rank, feature_idx in enumerate(sorted_indices):
                FI_indices[i, j, rank] = keys[feature_idx]

    # map FI indices back to grid
    fi_data = np.full((original_shape[0], original_shape[1], len(fnames)), -1, dtype=int)

    for i in range(num_domains_x):
        for j in range(num_domains_y):
            lat_slice = slice(i * domain_lat, (i + 1) * domain_lat)
            lon_slice = slice(j * domain_lon, (j + 1) * domain_lon)
            fi_data[lat_slice, lon_slice, :] = FI_indices[i, j, :]

    print("Feature importance data shape:", fi_data.shape)

    # aggregate feature importance scores over domains
    result_list = {}
    for i in range(num_domains_x):
        for j in range(num_domains_y):
            for key, val in fi_scores_array[i, j].items():
                result_list[key] = result_list.get(key, 0) + val
    print("Global feature importances:", result_list)

    # process R2 scores to original grid
    r2_map = np.full(original_shape, -1.0)
    for i in range(num_domains_x):
        for j in range(num_domains_y):
            lat_slice = slice(i * domain_lat, (i + 1) * domain_lat)
            lon_slice = slice(j * domain_lon, (j + 1) * domain_lon)
            r2_map[lat_slice, lon_slice] = r2_array[i, j]
    print("R2 scores map shape:", r2_map.shape)

    # process R2 scores to original grid
    rmse_map = np.full(original_shape, -1.0)

    # assign values to map
    for i in range(num_domains_x):
        for j in range(num_domains_y):
            lat_slice = slice(i * domain_lat, (i + 1) * domain_lat)
            lon_slice = slice(j * domain_lon, (j + 1) * domain_lon)
            rmse_map[lat_slice, lon_slice] = rmse_array[i, j]
    print("RMSE scores map shape:", rmse_map.shape)

    # reshape FI and R2 data for saving (FI data to (order,feat,lat,lon) and R2 and rmse to (value,lat,lon) )
    final_fi = np.expand_dims(np.transpose(fi_data, (2, 0, 1)), axis=0)
    r2_data = np.expand_dims(r2_map, axis=0)
    rmse_data = np.expand_dims(rmse_map, axis=0)

    print("final_fi shape:", final_fi.shape)

    log_tag = 'log_' if use_log else ''

    # file naming
    fi_filename = f"FI_{log_tag}{len(fnames)}feats_static_{domain_lat}_by_{domain_lon}_grids_w_{modelsel}.nc"
    r2_filename = f"R2s_for__{fi_filename}"
    rmse_filename = f"RMSE_for__{fi_filename}"

    # save results to netCDF
    write_4D_grid(os.path.join(output_dir, fi_filename), final_fi, lon, lat, range(len(fnames)), [1], 'FI index')
    write_2D_grid(os.path.join(output_dir, r2_filename), r2_data, lon, lat, 'R2 score')
    write_2D_grid(os.path.join(output_dir, rmse_filename), rmse_data, lon, lat, 'RMSE')

    # save results to png
    fi_filename_png = os.path.join(output_dir, fi_filename).replace('.nc','_map.png')
    r2_filename_png = os.path.join(output_dir, r2_filename).replace('.nc','_map.png')
    rmse_filename_png = os.path.join(output_dir, rmse_filename).replace('.nc','_map.png')

    rank = 1 
    for level in final_fi[0]:
        plot_map(level,lon,lat,fi_filename_png.replace('map',f'rank{rank}'),
            'Index of Feature','rainbow',fnames)
        rank+=1
    plot_map(r2_data[0],lon,lat,r2_filename_png, 'R²','cividis',fnames) # [0] for flattening the data
    plot_map(rmse_data[0],lon,lat,rmse_filename_png,'RMSE','viridis',fnames)

    print("Feature importance and R2/RMSE values saved!")

    # save metadata
    with open(os.path.join(output_dir, f'Lower_res_{modelsel}_metadata_{domain_lat}_by_{domain_lon}_grids.txt'), 'w') as f:    
        f.write(f"Mean R2 score: {np.mean(r2_data):.4f}\n")
        f.write(f"Mean RMSE score: {np.mean(rmse_data):.4f}\n")
        # dictionary with feature names and indices
        feature_names = {
            0: 'RH',
            1: 'MMR',
            2: 'Refrac',
            3: 'AbsCoef',
            4: 'CLT'
        }
        f.write(f"Summed feature importances (sorted by FI):\n")
        for key, importance in sorted(result_list.items(), key=lambda x: x[1], reverse=True):
            feature = feature_names.get(key, f"Feature_{key}")
            f.write(f" {feature} (index {key}): {importance}\n")

    end_time = time.time()-start_time
    print(f"FI study for lower resolution grid done in {end_time}")
    return

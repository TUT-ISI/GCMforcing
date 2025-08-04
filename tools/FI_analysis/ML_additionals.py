import pandas as pd
import numpy as np
from sklearn.metrics import r2_score, PredictionErrorDisplay, root_mean_squared_error
from sklearn.ensemble import RandomForestRegressor
import dask
from dask import delayed, compute
from dask.diagnostics import ProgressBar
from sklearn.inspection import permutation_importance
from scipy.stats import randint, uniform

import time
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import LinearRegression
# explicitly require this experimental feature
from sklearn.experimental import enable_halving_search_cv
# now you can import normally from model_selection
from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.ensemble import GradientBoostingRegressor

"""
by  Atte Laakso / Aalto University

Auxiliary functions for FI analysis
"""

# Function for adding coordinate indexes as features to data created in get_data workflow
def add_coordinates_to_features(lat, lon, x_train, x_test, fnames):
    num_lat = lat.size
    num_lon = lon.size
    grid_points_per_sample = num_lat * num_lon  # Total points in each grid
    
    # Create grid index arrays
    row_indices, col_indices = np.meshgrid(np.arange(num_lat), np.arange(num_lon), indexing='ij')
    
    # Flatten indices the same way data is flattened in get_data workflow
    row_flat = row_indices.flatten()
    col_flat = col_indices.flatten()
    
    # Normalize coordinates to interval of [0,1]
    row_flat_norm = row_flat / (num_lat - 1) if num_lat > 1 else row_flat
    col_flat_norm = col_flat / (num_lon - 1) if num_lon > 1 else col_flat

    # Calculate how many samples in train/test
    num_train_samples = len(x_train) // grid_points_per_sample
    num_test_samples = len(x_test) // grid_points_per_sample
    
    # Sanity check
    assert len(x_train) == num_train_samples * grid_points_per_sample, "Train data size mismatch"
    assert len(x_test) == num_test_samples * grid_points_per_sample, "Test data size mismatch"
    
    # Repeat coordinate indices for each sample
    row_coords_train = np.tile(row_flat_norm, num_train_samples)
    col_coords_train = np.tile(col_flat_norm, num_train_samples)
    row_coords_test = np.tile(row_flat_norm, num_test_samples)
    col_coords_test = np.tile(col_flat_norm, num_test_samples)
    
    # Convert inputs to DataFrames if not already
    x_train_df = pd.DataFrame(x_train, columns=fnames)
    x_test_df = pd.DataFrame(x_test, columns=fnames)
    
    # Add coordinate features
    x_train_df['row_idx'] = row_coords_train
    x_train_df['col_idx'] = col_coords_train
    x_test_df['row_idx'] = row_coords_test
    x_test_df['col_idx'] = col_coords_test
    
    # Extend feature name list
    fnames_extended = list(fnames) + ['row_idx', 'col_idx']
    
    return x_train_df, x_test_df, fnames_extended

# Function for dividing data to smaller domains
def form_gridwise_data(dataframe, original_shape, domain_lat, domain_lon):
    flattened_data = np.array(dataframe)
    print("dataset has shape of ", flattened_data.shape)
    
    if flattened_data.ndim < 2:
        raise ValueError("Input data must have more than one feature (column).")
    
    num_features = flattened_data.shape[1] # columns
    num_samples = flattened_data.shape[0] # rows

    # compute expected shape and checks
    total_elements = num_samples
    elements_per_sample = original_shape[0] * original_shape[1] # size of one map of data
    num_datasets = total_elements // elements_per_sample
    
    if total_elements % elements_per_sample != 0:
        raise ValueError("Total elements not divisible by original dataset shape")

    num_domains_lat = original_shape[0] // domain_lat
    num_domains_lon = original_shape[1] // domain_lon

    # initialize [i][j] structure
    domained_data = [[None for _ in range(num_domains_lon)] for _ in range(num_domains_lat)]

    for i in range(num_domains_lat):
        for j in range(num_domains_lon):
            domain_samples = []
            for n in range(num_datasets):
                # define start and end indices
                start = n * elements_per_sample
                end = (n + 1) * elements_per_sample
                # take sample out of all data
                sample = flattened_data[start:end].reshape(original_shape[0], original_shape[1], num_features)
                # restrict to domain
                domain = sample[i*domain_lat:(i+1)*domain_lat, j*domain_lon:(j+1)*domain_lon, :]
                domain_flat = domain.reshape(-1, num_features) # flatten the domain data
                domain_samples.append(domain_flat) # allocate domain data
            # stack samples along sample axis (here time axis)
            domained_data[i][j] = np.vstack(domain_samples)
    
    return domained_data

# Function for developing selected ML model for given data and for calculating 
# FI values and R2 and rmse scores
def get_FI(features, targets, test_feat, test_targ, fnames, rs, modelsel):

    assert features.shape[0] == targets.shape[0], "Mismatch in training samples"
    assert test_feat.shape[0] == test_targ.shape[0], "Mismatch in test samples"
    
    # ravel target datas if they have that attribute
    targets = targets.ravel() if hasattr(targets, 'ravel') else targets
    test_targ = test_targ.ravel() if hasattr(test_targ, 'ravel') else test_targ

    # model and hyperparameter sets (LR will be used if not one of these)
    model_configs = {
        "RF": (
            RandomForestRegressor(random_state=rs),
            {
                "n_estimators":     randint(100, 301),
                "max_depth":        [10, 20, None],
                "max_features":     ["sqrt", 0.5],
                "min_samples_split":randint(2, 6),
                "min_samples_leaf": randint(1, 5),
                "bootstrap":        [True],
                "criterion":        ["squared_error"]
            }
        ),
        "DT": (
            DecisionTreeRegressor(random_state=rs),
            {
                "max_depth":        [5, 10, 20, None],
                "min_samples_split":randint(2, 10),
                "min_samples_leaf": randint(1, 6),
                "max_features":     ["sqrt", None],
                "ccp_alpha":        uniform(0.0, 0.01)
            }
        ),
        "GB" : (
            GradientBoostingRegressor(random_state=rs),
            {
                "n_estimators": randint(200, 601),
                "learning_rate": uniform(0.01, 0.29),
                "max_depth": randint(2, 9),
                "min_samples_leaf": randint(1, 11),
                "subsample": uniform(0.6, 0.4)
            }
        )
    }

    # select model and it's config
    if modelsel in model_configs:
        base_model, param_dist = model_configs[modelsel]

        # quick hyperparameter tuning
        model_grid = HalvingRandomSearchCV(
            estimator=base_model,
            param_distributions=param_dist,
            cv=5, factor=3, n_jobs=1,
            random_state=rs,
            error_score='raise',
            min_resources=features.shape[0] // 2
        )
        model_grid.fit(features, targets)
        model = model_grid.best_estimator_
    else:
        # otherwise uses LR!
        model = LinearRegression()
        model.fit(features, targets)

    # feature importances with permutation (ensure that test samples are big enough)
    if test_feat.shape[0] > 10:
        result = permutation_importance(model, test_feat, test_targ, n_repeats=5, random_state=rs, n_jobs=1)
        fi_score = {i: result.importances_mean[i] for i in range(len(fnames))}
    else:
        fi_score = {i: 0.0 for i in range(len(fnames))}

    # calcualate also R2 and RMSE
    pred = model.predict(test_feat)
    r2 = r2_score(test_targ, pred)
    rmse = root_mean_squared_error(test_targ, pred)

    return fi_score, r2, rmse

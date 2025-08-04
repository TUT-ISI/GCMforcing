import os, sys
import math
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import random
from ML_eval import evaluate
from sklearn.ensemble import RandomForestRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import GradientBoostingRegressor
# explicitly require this experimental feature
from sklearn.experimental import enable_halving_search_cv
# now you can import normally from model_selection
from sklearn.model_selection import HalvingRandomSearchCV
import seaborn as sns
import joblib

def global_ml(x_train, y_train, x_test, y_test, original_shape, fnames, use_log, output_dir,
     rs, modelsel):
    
    """
    by Atte Laakso / Aalto University

    Develops defined ML model based on datasets given to it utilizing HalvingRandomSearchCV for
    parameter tuning. Calculates feature importances, R2 scores and rmse, then saves these results.

        Parameters:
        - x_train, y_train, x_test, y_test: pandas DataFrames or arrays
        - original_shape: tuple of (lat_size, lon_size) original grid dimensions
        - fnames: list of feature names
        - use_log: boolean, whether data is log transformed
        - rs: int, constant to feed to random states of Ml models
        - output_dir: string, directory for saving results
        - modelsel: string, defines model to be used (see code below + ML_basis)
        
        Returns: None but saves the calculated data to outputfiles using ML_eval.py
    """

    # keep track of time
    start_time = time.time()

    model_file = os.path.join(output_dir, f"trained_{modelsel}_model.pkl")

    if not os.path.exists(model_file):

        # Sample fraction of training data for tuning to lighten the program
        sample_frac = 0.2
        x_sample = x_train.sample(frac=sample_frac, random_state=rs)
        y_sample = y_train.loc[x_sample.index]  # match indices

        # verify dataset consistency
        assert len(x_test) == len(y_test), "Mismatch in test set lengths"
        assert len(x_sample) == len(y_sample), "Mismatch in test set lengths"
        print(f"Training set size: {len(x_train)} | Test set size: {len(x_test)}")
        print(f"Training sample size: {len(x_sample)}")

        print(f"Starting {modelsel} training with hyperparameter tuning...")
        ml_start_time = time.time()

        # initialize model and halving grid search
        if modelsel=='RF':
            model = RandomForestRegressor(random_state=rs)
            # define hyperparameter search space
            param_dist = {
                'n_estimators': [50, 100, 200],
                'max_depth': [10, 20, 30, None],
                'max_features': ['sqrt', 0.5],
                'min_samples_split': [2, 5, 10],
                'min_samples_leaf': [2, 4]
            }
        elif modelsel=='DT':
            model = DecisionTreeRegressor(random_state=rs)
            # define hyperparameter search space
            param_dist = {
                'max_depth': [5, 10, 20, None],
                'max_features': [None, 'sqrt', 0.7],
                'min_samples_split': [2, 5, 10],
                'min_samples_leaf': [1, 2, 4]
            }
        elif modelsel=='GB':
            model = GradientBoostingRegressor(random_state=rs)
            # define hyperparameter search space
            param_dist = {
                'learning_rate': [0.01, 0.05, 0.1, 0.2],
                'n_estimators': [100, 200, 500],
                'max_depth': [3, 5, 10],
                'max_features': [None,'sqrt', 0.7],
                'min_samples_split': [2, 5, 10],
                'min_samples_leaf': [1, 2, 4],
                'subsample': [0.6, 0.8, 1.0]
            }
        else:
            model = LinearRegression()
            # define hyperparameter search space
            param_dist = {}

        # use HalvingRandomSearchCV for finding optimal parameters for DT, GB and RF only
        if param_dist:
            model_grid = HalvingRandomSearchCV(
                estimator=model, param_distributions=param_dist,
                cv=3, factor=3, verbose=10, n_jobs=-1,
                random_state=rs, error_score='raise', min_resources=100000
            )
            model_grid.fit(x_sample, y_sample.values.ravel())
            best_model = model_grid.best_estimator_
            print("Best parameters found:")
            print(model_grid.best_params_)
        else:
            print("No hyperparameter tuning for model:", modelsel)
            best_model = model

        print("Fitting model with best parameters")
        # get best model on full data
        best_model.fit(x_train, y_train.values.ravel())

        ml_end_time = time.time()
        print(f"Training and tuning done. Took {ml_end_time - ml_start_time} seconds.")

        print("Saving the model")
        # save trained model
        joblib.dump(best_model, model_file)

    else:
        print("Reading model from file")
        best_model = joblib.load(model_file)

    # free some memory
    del x_train, y_train, original_shape
    
    # evaluate on test set, calculate R2, RMSE and FI values, plot them and save to output
    evaluate(best_model, x_test, y_test,output_dir,rs,modelsel)

    total_time = time.time() - start_time
    print(f"Total loc_ml execution time: {total_time} seconds.")
    return
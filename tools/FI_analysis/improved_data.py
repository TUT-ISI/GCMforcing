import os, sys
import math
import netCDF4 as nc
import numpy as np
import fnmatch
import pandas as pd
import random
from collections import defaultdict
from additionals import *
from feature_data import create_combos

"""
by Atte Laakso / UEF

Form ML compatible data structures from feature and label intended data stored in netCDF files

Parameters:
- models: string identifiers for the models to use
- path_to_goal_set: path to main directory that contains all label files (can include subdirs)
- fname_xxx: SALSA based data
- all_feats: list of features wanted to be included (eg. ['rh','mmr','ref_ind','abs_coef','clt'])

Returns:
- training and testing datasets for features and labels and length of data rows
"""


def create_data(models, path_to_goal_set, path_to_folder, fname_tracer, fname_vphysc, fname_ham, fname_activ, fname_echam, all_feats, use_ADRE):    
    
    x, y = pd.DataFrame(), pd.DataFrame(columns=['target'])
    x_test, y_test = pd.DataFrame(), pd.DataFrame(columns=['target'])
    batch_len = 0

    for model in models:
        print(f"\n--- Processing model: {model} ---")

        #---------------------------
        # Label (target) data file search
        #---------------------------
        a_or_d = 'ADRE' if use_ADRE else 'DELTA'
        goal_dictionary = defaultdict(list)
        target_paths = []

        for root, dirs, files in os.walk(path_to_goal_set):
            # Skip unwanted directories
            if any(skip in root for skip in ['old', 'temporary']):
                continue

            # Find matching files in current dir
            found_files = [
                os.path.join(root, f)
                for f in files
                if fnmatch.fnmatch(f, f"*{a_or_d}*") and os.path.isfile(os.path.join(root, f))
            ]

            if len(found_files) == 12:
                target_paths.append((root, found_files))
            
        print(f"Found {len(target_paths)} directories with valid target files")

        #---------------------------
        # Parse label data files into dictionary
        #---------------------------
        # Debugging
        total_files = 0
        skipped = 0
        duplicate_keys = 0
        
        for dir_path, files in target_paths:
            for file in files:
                total_files += 1
                
                # There are probably also other files in these dirs so do not include them
                if a_or_d not in file:
                    skipped +=1
                    continue
                
                # Search all feature keys you can find from the file name
                file_path = os.path.join(dir_path, file)
                found_keywords = sorted([kw for kw in all_feats if kw in file]) or ['_']

                # Extract [month]+[model name]+[features], if data from data_generation
                month = file.partition(a_or_d)[2][5:7]  # As it is formatted
                after_for = file.partition('for_')[2]  # Gives '<model_name>.nc'
                model_id = model_id = after_for[:-3]     # Strips off '.nc'
                key = f"{month}{model_id}?" + '!'.join(found_keywords) # Form key

                if key in goal_dictionary:
                    duplicate_keys += 1

                # Extract data from corresponding file
                with nc.Dataset(file_path, 'r') as ds:
                    var = ds.variables['ADRE'] if use_ADRE else ds.variables['dADRE']
                    data = np.array(var)

                # Store data under key
                goal_dictionary[key] = data

        print(f"Total files processed: {total_files}")
        print(f"Unique keys stored: {len(goal_dictionary)}")
        print(f"Duplicate keys encountered: {duplicate_keys}")
        print(f"Files skipped (missing '{a_or_d}'): {skipped}")

        #---------------------------
        # Combine features and labels
        #---------------------------
        comb = 1
        outcomes = [True, False]
        probabilities = [0.7, 0.3]

        for key, target_data in goal_dictionary.items():
            month = int(key[:2])
            model_key = key[2:].split('?')[0]

            if model_key not in model:
                continue
            print(f"\n Processing combination {comb} for model {model}")

            # Feature toggle based on key (feats after ? and separated from each other by !)
            raw_feats = key.split('?')[1].split('!')
            features_used = [f for f in raw_feats if f != '_']
            print("Features used:", features_used)

            param_combo = [
                any('rh' in f for f in features_used),
                any('mmr' in f for f in features_used),
                any('ref' in f for f in features_used),
                any('clt' in f for f in features_used)
            ]
            # NOTE: Extend this list if you add more features!

            # Get feature grids for this combination using function from feature_data.py
            try:
                feature_dict = create_combos(param_combo, model, path_to_folder, fname_ham, fname_echam, fname_vphysc, fname_tracer)
            except Exception as e:
                print(f"Failed to create features for key {key}: {e}")
                continue

            # Build feature arrays
            try:
                features_array = [feature_dict[feat][month-1].flatten() for feat in feature_dict.keys()]
                combined_features = np.vstack(features_array).T # Transpose to get data corresponding to each feature under columns
            except Exception as e:
                print(f"Error combining features for {key}: {e}")
                continue

            # Name columns and rows for pandas
            feature_rows = pd.DataFrame(combined_features, columns=feature_dict.keys())
            target_rows = pd.DataFrame(target_data[0].flatten(), columns=['target'])

            # Get the mask of valid (non-NaN) labels (targets)
            valid_mask = ~target_rows['target'].isna()

            # Create masked versions of the dataframes
            feature_rows_masked = feature_rows.copy()
            target_rows_masked = target_rows.copy()

            # Apply NaN to entire rows where the target is NaN to handle those later 
            # but to keep the shape of the data
            feature_rows_masked.loc[~valid_mask, :] = np.nan
            target_rows_masked.loc[~valid_mask, :] = np.nan

            # Count how many rows were masked
            removed_rows = (~valid_mask).sum()
            if removed_rows > 0:
                print(f"Masked {removed_rows} rows with NaNs in target")

            # Based on chances given in first part of this program, divide data to train and test
            if random.choices(outcomes, probabilities)[0]:
                x = pd.concat([x, feature_rows_masked], ignore_index=True)
                y = pd.concat([y, target_rows_masked], ignore_index=True)
            else:
                x_test = pd.concat([x_test, feature_rows_masked], ignore_index=True)
                y_test = pd.concat([y_test, target_rows_masked], ignore_index=True)
                batch_len = len(target_rows_masked)

            print(f"Combination {comb} complete | Features: {feature_rows_masked.shape} | Target: {target_rows_masked.shape}")
            comb += 1

    return x, y, x_test, y_test, batch_len


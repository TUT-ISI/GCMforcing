import os, sys
import numpy as np
import pandas as pd
from improved_data import create_data
from additionals import get_salsa_files
import time
current_directory = os.getcwd()
import logging
# initialize logging for later debugging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

"""
by Atte Laakso / UEF

Program calls method from improved_data.py to create features and labels and stores them
based on given parameters.

Run this before starting ML development in ML_basis as data needed there must be created here
"""
 
#---------------------------------------------------------------------------------------------------------------------------------------------------

# list of names for the other models in the comparison
models = ["CAM5-ATRAS","ECHAM6.3-HAM2.3","ECHAM6.3-SALSA2.0","GISS-ModelE2p1p1-OMA","GFDL-AM4"]

# define path to the dicitonary containing all the target sets
path_to_goal_set = '../output_DRE'

# define path to the AeroCom model dictionary
path_to_AeroCom = '../AeroCom_models'

# create list of all features (as they are stated in file names so eg. ['rh','mmr','refrac','clt'])
all_feats=['rh','mmr','refrac','clt']

# define path to SALSA data (fix)
sp = '../SALSA_2010'

# get salsa files
fname_tracer, fname_vphysc, fname_ham, fname_activ, fname_echam = get_salsa_files("vbs_sensitivity_base_2010.01" , sp)

# use direct radiative effects. (In other case the dADRE values will be used)
use_ADRE = False

# define datafiles prefixes (In-common part of the name of the files containing training and testing datasets)
ML_data = "dARE_ML_data_250721"


#---------------------------------------------------------------------------------------------------------------------------------------------------
# function for other programs to fetch the current prefix for data files
def prefix_name():
    return ML_data
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# main function

def main():
    # get salsa files
    fname_tracer, fname_vphysc, fname_ham, fname_activ, fname_echam = get_salsa_files("vbs_sensitivity_base_2010.01" , sp)

    # keep track of time
    start_time = time.time()

    # if data allready exists, read it from files (If you want new data, delete the files or change prefix)
    if (os.path.exists(f'{ML_data}_x_train.parquet') and
        os.path.exists(f'{ML_data}_y_train.parquet') and
        os.path.exists(f'{ML_data}_x_test.parquet') and
        os.path.exists(f'{ML_data}_y_test.parquet')):
        
        # data from files
        x_train = pd.read_parquet(f'{ML_data}_x_train.parquet')
        y_train = pd.read_parquet(f'{ML_data}_y_train.parquet')
        x_test = pd.read_parquet(f'{ML_data}_x_test.parquet')
        y_test = pd.read_parquet(f'{ML_data}_y_test.parquet')
        
        # other variables
        with open(f'{ML_data}_metadata.txt', 'r') as f:
            batch_len = int(f.readline().strip())

    else:
        # create data if some of the files was not found
        x_train, y_train, x_test, y_test, batch_len = create_data(models, path_to_goal_set, path_to_AeroCom, fname_tracer, fname_vphysc, fname_ham, fname_activ, fname_echam, all_feats, use_ADRE)
        
        # save DataFrames to parquet files
        x_train.to_parquet(f'{ML_data}_x_train.parquet')
        y_train.to_parquet(f'{ML_data}_y_train.parquet')
        x_test.to_parquet(f'{ML_data}_x_test.parquet')
        y_test.to_parquet(f'{ML_data}_y_test.parquet')
        
        # save other variables
        with open(f'{ML_data}_metadata.txt', 'w') as f:
            f.write(f"{batch_len}\n")
        
    print(f"Training sets are {len(x_train)}, {len(y_train)} and test sets are {len(x_test)} and {len(y_test)}")
    print(f"Training sets are shaped {x_train.shape}, {y_train.shape} and test sets are {x_test.shape} and {y_test.shape}")

    end_time = time.time()
    # print the duration
    print(f"initialization done. Initializing took {end_time - start_time} seconds")
    # +
    # get feature names
    fnames = x_train.columns
    print("Features are ",fnames)

if __name__ == "__main__":
    main()

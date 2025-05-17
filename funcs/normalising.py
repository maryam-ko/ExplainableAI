#!/bin/python

import pandas as pd
import numpy as np
import copy
from functools import reduce
import os


# ----------------- #

def create_dataframe_per_dataset(dataframe):
    # Ensure 'DatasetName' is set as the index
    if 'DatasetName' not in dataframe.index.names:
        dataframe.set_index('DatasetName', inplace=True)
    
    dataset_list = []
    for column in dataframe.columns:
        # Avoid creating unnecessary copies; work with the column directly
        df = dataframe[[column]].dropna()
        
        # Only append non-empty DataFrames
        if not df.empty:
            dataset_list.append(df)
        else:
            print(f"[INFO] Skipping empty column: {column}")
        
        # Free memory after processing each dataset to avoid holding large objects
        del df  # Delete the current dataframe to free up memory

    # Run garbage collection to ensure memory is released
    import gc
    gc.collect()
    
    return dataset_list

# ----------------- #

def MinMax_normalize_and_merge(dictionary, scalar): # MinMax normalise all dataframes in a dictionary and return a merged dataframe
    MinMax_dict = copy.deepcopy(dictionary) # copy dictionary
    for dataset in MinMax_dict:
        if 'DatasetName' not in MinMax_dict[dataset].index.names:
            MinMax_dict[dataset].set_index('DatasetName', inplace=True)
        MinMax_dict[dataset][MinMax_dict[dataset].columns] = scalar.fit_transform(MinMax_dict[dataset][MinMax_dict[dataset].columns]) # MinMax normalise each dataset
        MinMax_dict[dataset].reset_index(inplace=True) # reset index
    MinMax_itr = MinMax_dict.values() # convert dict values (dfs) to list
    MinMax_merged = reduce(lambda  left,right: pd.merge(left,right,on=['DatasetName'], how='outer'), MinMax_itr) # merge datasets on column
    return MinMax_merged

# ----------------- #

def create_and_save_per_dataset(dataframe, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    if 'DatasetName' not in dataframe.index.names:
        dataframe.set_index('DatasetName', inplace=True)

    for column in dataframe.columns:
        df = dataframe[[column]].dropna()
        if not df.empty:
            df.to_csv(os.path.join(output_dir, f"{column}.csv"))
        else:
            print(f"[INFO] Skipping empty column: {column}")

#!/bin/python

import pandas as pd
import numpy as np
import copy
from functools import reduce
import os


# ----------------- #

# def create_dataframe_per_dataset(dataframe): # create new dataframe for each column in a larger dataset
#     if 'DatasetName' not in dataframe.index.names:
#         dataframe.set_index('DatasetName', inplace=True)
        
#     dataset_list = [dataframe[[column]].dropna() for column in dataframe.columns]
#     return dataset_list

def create_dataframe_per_dataset(dataframe):
    if 'DatasetName' not in dataframe.index.names:
        dataframe.set_index('DatasetName', inplace=True)

    # Create one DataFrame per dataset (i.e., per row in the matrix)
    dataset_list = [
        group.dropna(axis=1)  # Drop phosphosites with NaN values for this dataset
        for _, group in dataframe.groupby(level='DatasetName')
    ]

    return dataset_list

# ----------------- #

# def MinMax_normalize_and_merge(dictionary, scalar): # MinMax normalise all dataframes in a dictionary and return a merged dataframe
#     MinMax_dict = copy.deepcopy(dictionary) # copy dictionary
#     for dataset in MinMax_dict:
#         if 'phosphosite_ID' not in MinMax_dict[dataset].index.names:
#             MinMax_dict[dataset].set_index('phosphosite_ID', inplace=True)
#         MinMax_dict[dataset][MinMax_dict[dataset].columns] = scalar.fit_transform(MinMax_dict[dataset][MinMax_dict[dataset].columns]) # MinMax normalise each dataset
#         MinMax_dict[dataset].reset_index(inplace=True) # reset index
#     MinMax_itr = MinMax_dict.values() # convert dict values (dfs) to list
#     MinMax_merged = reduce(lambda  left,right: pd.merge(left,right,on=['phosphosite_ID'], how='outer'), MinMax_itr) # merge datasets on column
#     return MinMax_merged

# def MinMax_normalize_and_merge(dictionary, scalar): # MinMax normalise all dataframes in a dictionary and return a merged dataframe
#     MinMax_dict = copy.deepcopy(dictionary) # copy dictionary
#     for dataset in MinMax_dict:
#         if 'DatasetName' not in MinMax_dict[dataset].index.names:
#             MinMax_dict[dataset].set_index('DatasetName', inplace=True)
#         if MinMax_dict[dataset].shape[0] > 0:
#             MinMax_dict[dataset][MinMax_dict[dataset].columns] = scalar.fit_transform(MinMax_dict[dataset][MinMax_dict[dataset].columns]) # MinMax normalise each dataset
#         else:
#             MinMax_dict[dataset][MinMax_dict[dataset].columns] = 0.0
#         MinMax_dict[dataset].reset_index(inplace=True) # reset index
#     MinMax_itr = MinMax_dict.values() # convert dict values (dfs) to list
#     MinMax_merged = reduce(lambda  left,right: pd.merge(left,right,on=['DatasetName'], how='outer'), MinMax_itr) # merge datasets on column
#     return MinMax_merged

def MinMax_normalize_and_merge(dictionary, scalar):
    MinMax_dict = {
        dataset: (
            df.set_index('DatasetName') if 'DatasetName' not in df.index.names else df
        ).pipe(lambda x: pd.DataFrame(
            scalar.fit_transform(x.select_dtypes(include='number')),  # only numeric
            index=x.index,
            columns=[f"{dataset}_{col}" for col in x.select_dtypes(include='number').columns]
        )).reset_index()
        for dataset, df in dictionary.items()
        if not df.select_dtypes(include='number').empty  # skip if no numeric data
    }

    if not MinMax_dict:
        raise ValueError("No numeric columns found in any dataset.")

    MinMax_merged = reduce(
        lambda left, right: pd.merge(left, right, on='DatasetName', how='outer'),
        MinMax_dict.values()
    )
    return MinMax_merged

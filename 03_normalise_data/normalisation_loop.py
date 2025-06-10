#!/bin/python

# ----------------- #
# LOAD DEPENDENCIES
# ----------------- #

import sys
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler

# ----------------- #
# IMPORT R-NORMALISED MATRIX 
# ----------------- #

matrix = pd.read_csv('/data/home/bt24990/ExplainableAI/03_normalise_data/MatrixCSVs/NBA-Matrix_Quantile.csv', header=0)
matrix = matrix.set_index('DatasetName')

# ----------------- #
# ENSURE OUTPUT DIRECTORY EXISTS
# ----------------- #
os.makedirs('./ScaledMatrix', exist_ok=True)

print(os.getcwd())

# ----------------- #
# START LOOPING THROUGH DATASETS
# ----------------- #
for dataset in matrix.index:
    # Get the row for the current dataset
    row = matrix.loc[dataset]
    
    # Convert the row to a numpy array and reshape it
    row_array = np.array(row).reshape(-1, 1)
    
    # Initialize the MinMaxScaler
    scaler = MinMaxScaler()
    
    # Fit and transform the data
    scaled_row = scaler.fit_transform(row_array)
    
    # Convert back to a DataFrame
    scaled_row_df = pd.DataFrame(scaled_row, index=row.index, columns=[dataset])
    
    # Save the scaled row to a CSV file
    safe_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in dataset)
    scaled_row_df.to_csv(f'./ScaledMatrix/{safe_name}_scaled.csv')


# ----------------- #
# READ BACK ALL SCALED FILES AND PRINT HEAD
# ----------------- #
for file in os.listdir('./ScaledMatrix/'):
    if file.endswith('_scaled.csv'):
        # Read the scaled CSV file
        df = pd.read_csv(f'./ScaledMatrix/{file}', index_col=0)
        
        # Print the first few rows of the DataFrame
        print(f"Contents of {file}:")
        print(df.head())
        print("\n")

# ----------------- #


all_dfs = []

for file in os.listdir('./ScaledMatrix/'):
    if file.endswith('_scaled.csv'):
        # Read the scaled CSV file
        df = pd.read_csv(f'./ScaledMatrix/{file}', index_col=0)

        # Count columns before filtering
        rows_before, cols_before = df.shape

        # Drop columns that are all NA
        df = df.loc[:, df.count() != 0]
        df = df.loc[:, (df != 0).any(axis=0)]

        # Count columns after filtering
        rows_after, cols_after = df.shape

        print(f"Original shape: {rows_before} rows × {cols_before} cols")
        print(f"Filtered shape: {rows_after} rows × {cols_after} cols")
        print(f"Removed: {rows_before - rows_after} rows, {cols_before - cols_after} columns")

        # Append to list
        all_dfs.append(df)

# Combine all DataFrames vertically (stack rows)
combined_df = pd.concat(all_dfs, axis=1)

combined_df  = combined_df.T   

combined_df = combined_df.reindex(sorted(combined_df.columns), axis=1)

# Save the combined DataFrame to CSV
combined_df.to_csv('/data/home/bt24990/ExplainableAI/03_normalise_data/MatrixCSVs/NormalisedMatrix-Zscore.csv', index=True, index_label="DatasetName")




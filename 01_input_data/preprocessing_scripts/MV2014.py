#!/bin/python

# ----------------- #
# LOAD DEPENDENCIES
# ----------------- #

import sys
import os
import pandas as pd
import numpy as np

grandparent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(grandparent_dir)

from funcs import preprocessing

# ----------------- #
# LOAD & CLEAN DATA
# ----------------- #

dataset = 'MV2014'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/MV2014_raw.csv', header=0)
print('Raw data loaded.')
print(f"Dataset Columns: {data.columns}")
print(data.head())  # Print first few rows to inspect data

# filter data to keep only those with localization probability >= 0.85
data = data[data['Localization prob'] >= 0.85] 

# Filtering out semi-colons from 'Amino acid', 'Positions within proteins', and 'Gene names' columns
data = data[~data['Amino acid'].str.contains(';', na=False)]
data = data[~data['Positions within proteins'].str.contains(';', na=False)]
data = data[~data['Gene names'].str.contains(';', na=False)]

# filter data
data['Sequence window'] = data['Sequence window'].str.replace('_', '')

print("Before renaming:", data.columns)  # Debugging
data = data.rename(columns={'Gene names': 'GeneName'})
print("After renaming:", data.columns)  # Verify if rename worked

preprocessing.match_seq_to_genename(data, 'Sequence window')
print('Amino acid sequences matched to gene names.')

data['Phosphosite'] = data['Amino acid'].astype(str) + '(' + data['Positions within proteins'].astype(str) + ')'

# Debugging: Check if 'Phosphosite' column is present after creation
if 'Phosphosite' not in data.columns:
    print("Error: 'Phosphosite' column not created!")
else:
    print("Phosphosite column created successfully.")

print(data[['Amino acid', 'Positions within proteins', 'Phosphosite']].head())  # Check values

# Keep only 'Phosphosite' and ratio columns
keepcols = ['Phosphosite'] + [col for col in data.columns if 'Ratio' in col]
data = data[keepcols]
print("Final dataset preview:")
print(data.head())  # Display first few rows

# log2 transform the ratio columns (Ratio columns)
ratio_columns = [col for col in data.columns if 'Ratio' in col]
data[ratio_columns] = data[ratio_columns].apply(pd.to_numeric, errors='coerce')

data = preprocessing.create_phos_ID(data) # call function to create phosphosite_ID column
print('Phosphosite IDs created.')

data = preprocessing.log2_transform(data)
print('Data has been log2 transformed.')


data = preprocessing.clean_phosID_col(data)

# Save only Phosphosite and the ratio columns to CSV
final_columns = ['Phosphosite'] + [col for col in data.columns if 'Ratio' in col]
data = data[final_columns]

data.to_csv(f'/Users/maryam/Documents/maryam-ko-QMUL-MSc-Project/PreprocessedDatasets/MV2014.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)
                              
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

preprocessing.match_seq_to_genename(data, 'Sequence window')
print('Amino acid sequences matched to gene names.')

print("Before creating Phosphosite, columns:", data.columns)

data['Phosphosite'] = data['Amino acid'].astype(str) + '(' + data['Positions within proteins'].astype(str) + ')'

# Debugging: Check if 'Phosphosite' column is present after creation
if 'Phosphosite' not in data.columns:
    print("Error: 'Phosphosite' column not created!")
else:
    print("Phosphosite column created successfully.")


print("After creating Phosphosite, columns:", data.columns)  # Check if 'Phosphosite' is added
print(data[['Amino acid', 'Positions within proteins', 'Phosphosite']].head())  # Check values

keepcols = ['Proteins', 'Gene names', 'Amino acid', 'Positions within proteins', 'Sequence window', 'Modified sequence', 
            'Localization prob', 'Phosphosite'] + [col for col in data.columns if 'Ratio' in col]
data = data[keepcols]

# log2 transform the ratios (Ratio columns)
ratio_columns = data.columns[2:12]  # Adjust this slice as needed based on your dataset
data[ratio_columns] = data[ratio_columns].astype(float)

data.rename(columns={'Gene names': 'GeneName'}, inplace=True)

data = preprocessing.create_phos_ID(data) # call function to create phosphosite_ID column
print('Phosphosite IDs created.')

data = preprocessing.log2_transform(data)
print('Data has been log2 transformed.')


data = preprocessing.clean_phosID_col(data)

data.to_csv(f'/Users/maryam/Documents/maryam-ko-QMUL-MSc-Project/PreprocessedDatasets/MV2014.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)
                              
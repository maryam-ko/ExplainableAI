#!/bin/python

# ----------------- #
# LOAD DEPENDENCIES
# ----------------- #

import sys
import os
import pandas as pd
import numpy as np

grandparent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(grandparent_dir)
from funcs import preprocessing 

# ----------------- #
# LOAD & CLEAN DATA
# ----------------- #

dataset = 'MV2014'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/Users/maryam/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/MV2014_raw.csv', header=0)
print('Raw data loaded.')

# filter data to keep only those with localization probability >= 0.85
data = data[data['Localization prob'] >= 0.85] 

# filter data
data['Sequence window'] = data['Sequence window'].str.replace('_', '')

preprocessing.match_seq_to_genename(data, 'Sequence window')
print('Amino acid sequences matched to gene names.')

data['Phosphosite'] = data['Gene names'] + "_" + data['Amino acid'] + "(" + data['Positions within proteins'].astype(str) + ")"

# filter data to only include rows with a single protein (no semicolons in columns)
data = data[
    ~data['Proteins'].str.contains(';', na=False) & 
    ~data['Gene names'].str.contains(';', na=False) & 
    ~data['Positions within proteins'].str.contains(';', na=False)
]
print(f'Rows with semicolons removed. Shape: {data.shape}')


# keep only the necessary columns 
keepcols = [
    'Proteins', 'Gene names', 'Amino acid', 'Positions within proteins', 'Sequence window', 'Modified sequence', 
    'Localization prob',
    'Ratio M/L normalized A1___1 (2 min)', 'Ratio M/L normalized A1___2 (2 min)', 'Ratio M/L normalized A1___3 (2 min)',
    'Ratio H/M normalized A2___1 (2 min)', 'Ratio H/M normalized A2___2 (2 min)', 'Ratio H/M normalized A2___3 (2 min)',
    'Ratio L/H normalized A3___1 (2 min)', 'Ratio L/H normalized A3___2 (2 min)', 'Ratio L/H normalized A3___3 (2 min)',
    'Ratio H/L normalized A1___1 (10 min)', 'Ratio H/L normalized A1___2 (10 min)', 'Ratio H/L normalized A1___3 (10 min)',
    'Ratio L/M normalized A2___1 (10 min)', 'Ratio L/M normalized A2___2 (10 min)', 'Ratio L/M normalized A2___3 (10 min)',
    'Ratio M/H normalized A3___1 (10 min)', 'Ratio M/H normalized A3___2 (10 min)', 'Ratio M/H normalized A3___3 (10 min)',
    'Ratio M/L normalized B1___1 (5 min)', 'Ratio M/L normalized B1___2 (5 min)', 'Ratio M/L normalized B1___3 (5 min)',
    'Ratio H/M normalized B2___1 (5 min)', 'Ratio H/M normalized B2___2 (5 min)', 'Ratio H/M normalized B2___3 (5 min)',
    'Ratio L/H normalized B3___1 (5 min)', 'Ratio L/H normalized B3___2 (5 min)', 'Ratio L/H normalized B3___3 (5 min)',
    'Ratio H/L normalized B1___1 (30 min)', 'Ratio H/L normalized B1___2 (30 min)', 'Ratio H/L normalized B1___3 (30 min)',
    'Ratio L/M normalized B2___1 (30 min)', 'Ratio L/M normalized B2___2 (30 min)', 'Ratio L/M normalized B2___3 (30 min)',
    'Ratio M/H normalized B3___1 (30 min)', 'Ratio M/H normalized B3___2 (30 min)', 'Ratio M/H normalized B3___3 (30 min)'
]
data = data[keepcols]
print(f'Columns filtered. Shape: {data.shape}')

# log2 transform the ratios (Ratio columns)
data.iloc[:, 7:] = data.iloc[:, 7:].astype(float)  # convert from string to float


# create phosphosite ID 
data = preprocessing.create_phos_ID(data)  
print('Phosphosite IDs created.')

data = preprocessing.log2_transform(data)  
print('Data has been log2 transformed.')

data = preprocessing.clean_phosID_col(data)

output_path = '/Users/maryam/Documents/maryam-ko-QMUL-MSc-Project/preprocessed_datasets/MV2014.csv'

data.to_csv(output_path, index=False)

print(f'{dataset} processed data has been saved to CSV successfully at {output_path}')


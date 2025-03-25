import sys
import os
import pandas as pd
import numpy as np

grandparent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(grandparent_dir)

from funcs import preprocessing

dataset = 'KS2014'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/mmc3.csv', header=0)
print('Raw data loaded.')

print(f"Dataset Columns: {data.columns}")
print(data.head())  # Print first few rows to inspect data

# filter data to keep only those with localization probability >= 0.85
data = data[data['Localization prob'] >= 0.85] 

# Filtering out semi-colons from 'Amino acid', 'Positions', and 'Gene names' columns
data = data[~data['Amino acid'].str.contains(';', na=False)]
data = data[~data['Gene names'].str.contains(';', na=False)]

# filter data
data['Sequence window'] = data['Sequence window'].str.replace('_', '')

preprocessing.match_seq_to_genename(data, 'Sequence window')
print('Amino acid sequences matched to gene names.')

data['Phosphosite'] = data['Amino acid'].astype(str) + '(' + data['Positions within proteins'].astype(str) + ')'

# Keep only 'Phosphosite' and ratio columns
keepcols = ['Phosphosite'] + ['GeneName'] + [col for col in data.columns if 'Intensity' in col]
data = data[keepcols]

print("Data after subsetting columns:", data)
print("Cols after subsetting:", data.columns)

# log2 transform the ratio columns 
Intensity_columns = [col for col in data.columns if 'Intensity' in col]
data[Intensity_columns] = data[Intensity_columns].apply(pd.to_numeric, errors="coerce")
data[Intensity_columns] = np.log2(data[Intensity_columns] + 1)  # Avoid log(0) errors
print("After transformation:")
print(data.head())  # Show the first few rows after processing

data = preprocessing.create_phos_ID(data) # call function to create phosphosite_ID column
print('Phosphosite IDs created.')

data = preprocessing.clean_phosID_col(data)
print("After cleaning phosphosite_ID column:")
print(data.head())

final_columns = ['phosphosite_ID'] + [col for col in data.columns if 'Intensity' in col]
data = data[final_columns]
print("Final dataset preview:")
print(data.head())  # Display first few rows
print(data.tail())  # Display last few rows

data.to_csv(f'/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/PreprocessedDatasets/KS2014.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)
                              
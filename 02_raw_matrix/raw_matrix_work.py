#!/bin/python

# ----------------- #
# LOAD DEPENDENCIES
# ----------------- #
print('About to import dependencies')
import sys
print('Sys loaded')
import os
print('OS loaded')
import pandas as pd
print('Pandas loaded')
import numpy as np
print('Numpy loaded')

grandparent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(grandparent_dir)
from funcs import preprocessing


# ----------------- #
# GENERATE RAW MATRIX HEADER
# ----------------- #

# load preprocessed datasets
# stores names of processed datasets
file_names = ['AST2020']

print('File names loaded')

files_dict = preprocessing.create_dict_per_dataset(file_names)

matrix_cols = preprocessing.create_matrix_header(files_dict)
print(f'Matrix header:', matrix_cols)



# ----------------- #
# LOAD MATRIX HEADER
# ----------------- #

matrix = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/02_raw_matrix/RawMatrixProcessing/raw-matrix-header.csv', header = 0)



# ----------------- #
# LOAD DATASET NAMES
# ----------------- #

# define names for each dataset
AST2020_names = ('AST2020_Control',	'AST2020_Control.1', 'AST2020_Control.2',	'AST2020_anti-CD3',	'AST2020_anti-CD3.1',	'AST2020_anti-CD3.2',	'AST2020_anti-CD3+PDL2',	'AST2020_anti-CD3+PDL2.1',	'AST2020_anti-CD3+PDL2.2')

 
# ----------------- #
# PAIR DATASET NAMES WITH FILENAMES
# ----------------- #

files_datasets = [('AST2020', AST2020_names) ]

intermed_matrix = preprocessing.add_rows_to_matrix(matrix, files_datasets, )
print(intermed_matrix)



# ----------------- #
# FORMAT MATRIX
# ----------------- #

# reorder matrix columns
cols = intermed_matrix.columns.tolist()
cols = cols[-2:-1] + cols[:-2]
raw_matrix = intermed_matrix[cols]

# convert columns to numeric
numeric_cols = [i for i in raw_matrix.columns if i not in ['DatasetName']]
for col in numeric_cols:
    raw_matrix.loc[:, col]=pd.to_numeric(raw_matrix[col])
    
# remove infinity values
raw_matrix = raw_matrix.replace([np.inf, -np.inf], np.nan)

# save raw matrix
raw_matrix.to_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/RawMatrix.csv', index=False)
print(f'Raw matrix saved successfully!', raw_matrix)


# ----------------- #

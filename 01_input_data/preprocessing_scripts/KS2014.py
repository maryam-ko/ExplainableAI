import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

dataset = 'KS2014'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/mmc3.csv', header=0)
print('Raw data loaded.')
data

data.columns = data.columns.str.strip()

# filter data to keep only those with localization probability >= 0.85
data = data[data['Localization prob'] >= 0.85]

# Filtering out semi-colons from 'Amino acid', 'Positions', and 'Gene names' columns
data = data[~data['Amino acid'].str.contains(';', na=False)]
data = data[~data['Gene names'].str.contains(';', na=False)]


# filter data
data['Sequence window'] = data['Sequence window'].str.replace('_', '')

def match_seq_to_genename(dataset, seq_column):
    '''
    Maps amino acid sequences to gene names using the loaded fasta file.
    
    args:
    =====
    dataset: <pd.Dataframe> with a column of amino acid sequences
    seq_column: <str> column name containing amino acid sequences
    
    out:
    ====
    dataset: <pd.Dataframe> with an additional column containing gene names
    '''    
    
    fasta_sequence = list(SeqIO.parse(f"/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta", "fasta"))
    
    gene_dict = {}
    
    # iterate over rows in seq_column
    for i in dataset[seq_column]:
        print(i)
        i_str = str(i)
        for seq_record in fasta_sequence:
            matches = re.findall(i_str, str(seq_record.seq))
            if matches:
                print(f"Match found for sequence: {seq_record}")
                gene_name_match = re.search(r"GN=(\w+)", seq_record.description)
                print('Gene name match:', gene_name_match)
                if gene_name_match:
                    gene_name = gene_name_match.group(1)
                    gene_dict[i] = gene_name
                    print(f"Match found: {i_str} -> {gene_name}")
                else:
                    print(f"No gene name match found in description for sequence: {i_str}")

    # map sequences to gene names           
    dataset['GeneName'] = dataset[seq_column].map(gene_dict) 
    print('Amino acid sequences matched to gene names.')

    return dataset

data = match_seq_to_genename(data, 'Sequence window')

data

# Ensure GeneName exists before proceeding
if 'GeneName' not in data.columns:
    raise ValueError("GeneName column is missing! Check match_seq_to_genename function.")

data['Phosphosite'] = data['Amino acid'].astype(str) + '(' + data['Position'].astype(str) + ')'

# Debugging: Check if 'Phosphosite' column is present after creation
if 'Phosphosite' not in data.columns:
    print("Error: 'Phosphosite' column not created!")
else:
    print("Phosphosite column created successfully.")

data

# Keep only 'Phosphosite' and ratio columns
keepcols = ['Phosphosite'] + ['GeneName'] + [col for col in data.columns if 'Intensity' in col]
data = data[keepcols]
data

print("Data after subsetting columns:", data)

print("Cols after subsetting:", data.columns)

# log2 transform the ratio columns 
Intensity_columns = [col for col in data.columns if 'Intensity' in col]

data[Intensity_columns] = data[Intensity_columns].apply(pd.to_numeric, errors="coerce")

def log2_transform(dataset):
    '''
    Log2 transform a dataset.
    
    args:
    =====
    dataset: <pd.Dataframe>
    
    out:
    ====
    dataset: <pd.Dataframe> with log2 transformed values

    '''
    cols_to_transform = dataset[Intensity_columns]
    dataset[Intensity_columns] = cols_to_transform.apply(np.log2)
    print('Data has been log2 transformed.')
    return dataset

data = log2_transform(data)
print(f"DataFrame after log2 transformation:\n{data}") # Print the DataFrame after log2 transformation


def create_phos_ID(dataset):
    '''
    Concatenates GeneName and Phosphosite columns.
    
    args:
    =====
    dataset: <pd.Dataframe> with columns 'GeneName' and 'Phosphosite'
    
    out:
    ====
    dataset: <pd.Dataframe> with 'phosphosite_ID' column and 'GeneName' + 'Phosphosite' columns dropped
    '''
    dataset.loc[:, 'phosphosite_ID'] = dataset['GeneName'].astype(str) + '_' + dataset['Phosphosite'].astype(str)
    dataset = dataset.drop(columns=['Phosphosite', 'GeneName'])
    print('Phosphosite IDs created.')
    return dataset

data = create_phos_ID(data) # call function to create phosphosite_ID column

print('Phosphosite IDs created.')

def clean_phosID_col(data):
    data = data[~data.phosphosite_ID.str.contains('nan', case = False)]
    data = data[~data.phosphosite_ID.str.contains(';', case = False)] # remove rows containing ';' in phosphosite_ID column
    data = data[~data.phosphosite_ID.str.contains('-', case = False)] # remove rows containing '-' in phosphosite_ID column
    
    # check whether there are any phosphosites with multiple measurements
    data_grouped = data.groupby(by = 'phosphosite_ID')
    if len(data) != len(data_grouped):
        data = data_grouped.mean()
        data.reset_index(inplace=True) # reset index
        print('Phosphosites with multiple measurements have been averaged')
    else:
        print('There are no phosphosites with multiple measurements')
        
    print(data)
        
    data = data.replace([np.inf, -np.inf], np.nan)
        
    if data.columns[0] != 'phosphosite_ID':
        phosphosite_ID = data.pop('phosphosite_ID')
        data.insert(0, 'phosphosite_ID', phosphosite_ID)
    return data

data = clean_phosID_col(data)
print("After cleaning phosphosite_ID column:")
data

data.to_csv(f'/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/PreprocessedDatasets/KS2014.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)
                              
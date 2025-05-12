import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

dataset = 'CK2022'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/pmic13537-sup-0010-suppmats8.csv', header=0)
print('Raw data loaded.')
data

data.columns = data.columns.str.strip()

# filter data to keep only those with localization probability >= 0.85
data = data[data['Localization probability'] >= 0.85] 

# Filtering out semi-colons from 'Amino acid', 'Positions within proteins', and 'Gene names' columns
data = data[~data['Amino acid'].str.contains(';', na=False)]
data = data[~data['Position'].astype(str).str.contains(';', na=False)]
data = data[~data['Gene names'].str.contains(';', na=False)]
data

# filter data
data['Sequence window'] = data['Sequence window'].str.replace('_', '')
data

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

    fasta_sequence = list(SeqIO.parse(open(f'/data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta'), "fasta"))
    
    
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
                # gene_name_match = re.search("GN=(\w+)", seq_record.description)
                if gene_name_match:
                    gene_name = gene_name_match.group(1)
                    gene_dict[i] = gene_name
                    print(f"Match found: {i_str} -> {gene_name}")
                else: 
                    print(f"No gene name found in description for sequence: {i_str}")
    
    # map sequences to gene names           
    dataset['GeneName'] = dataset[seq_column].map(gene_dict) 
    print('Amino acid sequences matched to gene names.')
    return dataset 
    
data = match_seq_to_genename(data, 'Sequence window')

data['Phosphosite'] = data['Amino acid'].astype(str) + '(' + data['Position'].astype(str) + ')'

keepcols = ['Phosphosite'] + ['GeneName'] + [col for col in data.columns if '-' in str(col)]
data = data[keepcols]
data

Ratio_columns = [col for col in data.columns if '-' in str(col)]
data[Ratio_columns] = data[Ratio_columns].apply(pd.to_numeric, errors='coerce')

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
    cols_to_transform = dataset[Ratio_columns]
    dataset[Ratio_columns] = cols_to_transform.apply(np.log2)
    print('Data has been log2 transformed.')
    return dataset

data = log2_transform(data)
print(f"DataFrame after log2 transformation:\n{data}") # Print the DataFrame after log2 transformation
data

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
    data = data[~data.phosphosite_ID.str.contains('nan', case=False, na=False)]
    data = data[~data.phosphosite_ID.str.contains(';', case=False, na=False)]
    data = data[~data.phosphosite_ID.str.contains('-', case=False, na=False)]

    # Add this line to remove decimals from phosphosite_ID (e.g., S123.0 -> S123)
    data['phosphosite_ID'] = data['phosphosite_ID'].apply(lambda x: re.sub(r'\((\d+)\.0+\)', r'(\1)', x))
    
    data_grouped = data.groupby(by='phosphosite_ID')
    
    if len(data) != len(data_grouped):
        numeric_cols = data.select_dtypes(include=[np.number]).columns.tolist()
        non_numeric_cols = data.columns.difference(numeric_cols + ['phosphosite_ID']).tolist()
        data_numeric = data_grouped[numeric_cols].mean()
        data_categorical = data_grouped[non_numeric_cols].first().reset_index()
        
        # Merge numeric and non-numeric parts
        data = pd.merge(data_categorical, data_numeric, on='phosphosite_ID')
        print('Phosphosites with multiple measurements have been averaged')
    else:
        print('There are no phosphosites with multiple measurements')

    # Replace inf values with NaNs
    data = data.replace([np.inf, -np.inf], np.nan)
    
    # Ensure phosphosite_ID is first column
    if data.columns[0] != 'phosphosite_ID':
        phosphosite_ID = data.pop('phosphosite_ID')
        data.insert(0, 'phosphosite_ID', phosphosite_ID)

    return data

data = clean_phosID_col(data)
print("After cleaning phosphosite_ID column:")
data

data.to_csv(f'/data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/PreprocessedDatasets/CK2022.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)


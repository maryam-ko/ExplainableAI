import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

dataset = 'LATM2013'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/c3mb25524g_7.csv', header=0)
print('Raw data loaded.')
data

data.columns = data.columns.str.strip()

data['Amino acid'] = data['Sites \nAcc. # 1'].str.extract(r'([A-Z])')        # Amino acid
data['Position'] = data['Sites \nAcc. # 1'].str.extract(r'(\d+)')    # Position
print(data[['Sites \nAcc. # 1', 'Amino acid', 'Position']].head())

data = data.rename(columns={'Gene Names\nAcc. # 1': 'Gene Names'})

# Filtering out semi-colons from 'Amino acid', 'Positions within proteins', and 'Gene names' columns
data = data[~data['Amino acid'].str.contains(';', na=False)]
data = data[~data['Position'].str.contains(';', na=False)]
data = data[~data['Gene Names'].str.contains(';', na=False)]
data = data[data['Gene Names'].str.strip().str.count(r'\s') == 0]
data = data[~data['Gene Names'].str.strip().eq('Deleted')]
data

# filter data
data['Phosposite Sequence Context'] = data['Phosposite Sequence Context'].str.replace('_', '')
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

    fasta_sequence = list(SeqIO.parse(open(f'/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta'), "fasta"))
    
    
    gene_dict = {}
    
    # iterate over rows in seq_column
    for i in dataset[seq_column]:
        print(i)
        i_str = str(i)
        for seq_record in fasta_sequence:
            matches = re.findall(re.escape(i_str), str(seq_record.seq))
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

data = match_seq_to_genename(data, 'Phosposite Sequence Context')

data['Phosphosite'] = data['Amino acid'].astype(str) + '(' + data['Position'].astype(str) + ')'

keepcols = ['Phosphosite'] + ['GeneName'] + [col for col in data.columns if 'Log2 Ratio' in col]
data = data[keepcols]
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

data.to_csv(f'/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/PreprocessedDatasets/LATM2013.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)



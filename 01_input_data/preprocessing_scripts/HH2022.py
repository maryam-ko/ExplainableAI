import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

dataset = 'HH2022'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/41598_2022_8430_MOESM12_ESM.csv', header=0)
print('Raw data loaded.')
data

filtered_data = data[data['Modifications'].str.contains(r'^Phospho \(STY\)$', case=False, na=False)]
data = filtered_data[filtered_data['Modifications'].str.count(',') == 0]

print(data[['Gene Names', 'Modifications', 'Modified sequence']])

def has_high_phospho_probability(prob_str):
    # Extract all numbers within parentheses, 
    probabilities = re.findall(r'\((\d+\.\d+)\)', prob_str)
    
    # Convert strings to float and check if any value is greater than 0.85
    return any(float(prob) > 0.85 for prob in probabilities)

data = filtered_data[filtered_data['Phospho (STY) Probabilities'].apply(has_high_phospho_probability)]

print(data[['Gene Names', 'Modified sequence', 'Phospho (STY) Probabilities']])


def extract_modified_amino_acid(modified_seq):
    modified_amino_acids = re.findall(r'([A-Z])\(ph\)', modified_seq)  # Only match phosphorylation
    return [aa for aa in modified_amino_acids]


data['Amino Acid'] = data['Modified sequence'].apply(extract_modified_amino_acid)

print(data[['Gene Names', 'Modified sequence', 'Phospho (STY) Probabilities', 'Amino Acid']])

def find_position_in_gene(dataset, seq_column, fasta_path):
    positions_dict = {}
    
    # Parse the FASTA file
    fasta_sequence = list(SeqIO.parse(open(fasta_path), "fasta"))
    
    # Iterate over rows in the dataset
    for i, row in dataset.iterrows():
        gene_sequence = row[seq_column]  # Get the sequence from the dataset
        
        for seq_record in fasta_sequence:
            # Check if the sequence is part of the gene in the FASTA record
            start_position = str(seq_record.seq).find(gene_sequence)
            
            # If a match is found, store the position
            if start_position != -1:
                positions_dict[row['Gene Names']] = start_position
                break  # Stop once the position is found for the first matching gene
    
    # Add the start positions to the dataset as a new column
    dataset['StartPosition'] = dataset['Gene Names'].map(positions_dict)
    return dataset

fasta_path = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta'
data = find_position_in_gene(data, 'Sequence', fasta_path)

# Print the data with start positions
print(data[['Gene Names', 'Sequence', 'StartPosition']])

# Select the necessary columns: Gene Names, Modified Amino Acids, Sequence, and Reporter Intensity columns
columns_to_display = ['Gene Names', 'Amino Acid', 'StartPosition', 'Sequence'] + [col for col in data.columns if 'Reporter intensity' in col]

# Display the data with the selected columns
print(data[columns_to_display])

# Filtering out semi-colons from 'Amino acid', 'Positions within proteins', and 'Gene names' columns
data = data[~data['Amino Acid'].str.contains(';', na=False)]
data = data[~data['Gene Names'].str.contains(';', na=False)]
data

# filter data
data['Sequence'] = data['Sequence'].str.replace('_', '')
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
    
data = match_seq_to_genename(data, 'Sequence')

data['Phosphosite'] = data['Amino Acids'].astype(str) + '(' + data['StartPosition'].astype(str) + ')'

keepcols = ['Phosphosite'] + ['GeneName'] + [col for col in data.columns if 'Reporter intensity' in col]
data = data[keepcols]
data

Intensity_columns = [col for col in data.columns if 'Reporter intensity' in col]
data[Intensity_columns] = data[Intensity_columns].apply(pd.to_numeric, errors='coerce')

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

data.to_csv(f'/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/PreprocessedDatasets/HH2022.csv', index=False)


print(dataset, 'has been saved to CSV successfully!', data)


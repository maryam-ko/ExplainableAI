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

def match_seq_to_genename(dataset, seq_column, chunk_size=1000):
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

    # Break the dataset into chunks to avoid overloading memory
    for chunk_start in range(0, len(dataset), chunk_size):
        chunk_end = min(chunk_start + chunk_size, len(dataset))
        chunk = dataset.iloc[chunk_start:chunk_end]
    
    # iterate over rows in seq_column
    for i in chunk[seq_column]:
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

        # After processing each chunk, map sequences to gene names for that chunk
        dataset.loc[chunk_start:chunk_end, 'GeneName'] = chunk[seq_column].map(gene_dict)
    

    # map sequences to gene names           
    dataset['GeneName'] = dataset[seq_column].map(gene_dict) 
    print('Amino acid sequences matched to gene names.')
    return dataset


data = match_seq_to_genename(data, 'Sequence window')

data
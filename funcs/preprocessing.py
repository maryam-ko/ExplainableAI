#!/bin/python

import pandas as pd
import numpy as np
import os
from Bio import SeqIO
import re
from functools import reduce


# ----------------- #

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
    
    fasta_file = '/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta'

    # Check if the file exists
    if not os.path.exists(fasta_file):
        print(f"File not found: {fasta_file}")
        return dataset
    
    print(f"Opening fasta file: {fasta_file}")
    fasta_sequence = list(SeqIO.parse(open(fasta_file), "fasta"))
    
    gene_dict = {}
    
    # iterate over rows in seq_column
    for i in dataset[seq_column]:
        i_str = str(i)
        for seq_record in fasta_sequence:
            matches = re.findall(i_str, str(seq_record.seq))
            if matches:
                gene_name_match = seq_record.description.split(' ')[1].split(' ')[0]
                # gene_name_match = re.search("GN=(\w+)", seq_record.description)
                if gene_name_match:
                    gene_dict[i] = gene_name_match
    
    # map sequences to gene names           
    dataset['GeneName'] = dataset[seq_column].map(gene_dict) 
    print('Amino acid sequences matched to gene names.')
    
    return dataset


# ----------------- #

def find_position_in_gene(dataset, seq_column):
    positions_dict = {}
    
    fasta_sequence = list(SeqIO.parse(open('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta'), "fasta"))

    # iterate over rows in the Sequence Window column of GG2009
    for i in dataset['Sequence']:
        # iterate over entries (gene ids and sequences) in the fasta_sequence object
        for seq_record in fasta_sequence:
            # find matches between the sequence in the dataset and the sequence in the fasta file
            matches = re.findall(i, str(seq_record.seq))
            # if matches are found, print the gene id, the sequence, and the length of the sequence
            if matches:
                # find the position of i within the seq_record.seq
                start_position = str(seq_record.seq).find(i)
                #print("Starting position of i within seq_record.seq:", start_position)
                positions_dict[i] = start_position # add start position to dictionary

    dataset['StartPosition'] = dataset[seq_column].map(positions_dict)
    
    
    
# ----------------- #

def get_position_and_gene(dataset, seq_column, position_column):
    gene_dict = {}
    residues_dict = {}

    fasta_sequence = list(SeqIO.parse(open('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/UP000005640_9606.fasta'), "fasta"))
        
    # get the gene name and amino acid from the fasta file
    for index, row in dataset.iterrows():  # iterate over rows in the DataFrame
        sequence = row[seq_column]
        position = row[position_column] - 1  # Python uses 0-based indexing
        for seq_record in fasta_sequence:
            if sequence in str(seq_record.seq):  # if sequence matches
                gene_name_match = seq_record.description.split(' ')[1].split(' ')[0]
                # gene_name_match = re.search("GN=(\w+)", seq_record.description)
                if gene_name_match:
                    gene_dict[sequence] = gene_name_match
                    # Translate the sequence into amino acids
                    protein_sequence = seq_record.seq
                    # Find the amino acid at the given position
                    residue = protein_sequence[position]
                    residues_dict[sequence] = residue

    dataset['GeneName'] = dataset['Sequence'].map(gene_dict)  # map sequences to gene names
    dataset['Residue'] = dataset['Sequence'].map(residues_dict)  # map sequences to amino acids
    
    

# ----------------- #
    
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
    print('Phosphosite IDs created.')
    return dataset



# ----------------- #

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
    cols_to_transform = dataset.columns.drop('phosphosite_ID')
    
    dataset[cols_to_transform] = dataset[cols_to_transform].astype(float).apply(np.log2)
    
    print('Data has been log2 transformed.')
    return dataset.set_index('phosphosite_ID')  # Ensure phosphosites are the index


# ----------------- #

def rename_col_by_index(dataframe, index_to_name_mapping):
    '''
    Rename columns by index
    
    args:
    =====
    dataframe: <pd.Dataframe>
    index_to_name_mapping: <dict> mapping of index (key) to column name (str value)
    
    out:
    ====
    dataframe: <pd.Dataframe> with renamed columns
    '''
    dataframe.columns = [index_to_name_mapping.get(i, col) for i, col in enumerate(dataframe.columns)]
    return dataframe



# ----------------- #

def get_ens_dict(file_path): 
    '''
    Create dictionary of gene names and gene IDs from Ensembl gtf file
    '''
    with open(file_path) as f:
        gtf = [x for x in f if not x.startswith('#') and 'gene_id "' in x and 'gene_name "' in x]
    if len(gtf) == 0:
        print('you need to change gene_id ' and 'gene_name "')
    gtf = dict(set(map(lambda x: (x.split('gene_id "')[1].split('"')[0], x.split('gene_name "')[1].split('"')[0]), gtf)))
    print('number of ID pairs retrieved from database: ', len(gtf))
    return gtf

    # explanation of function at: https://www.youtube.com/watch?v=ve_BoDw1s7I



# ----------------- #

def create_dict_per_dataset(file_names):
    files_dict = {}
    for file in file_names:
        files_dict[file] = pd.read_csv(f'/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/data_files/PreprocessedDatasets/{file}.csv', header=0)
        print(f"{file} added to dict")
    print('Datasets have been loaded into dictionary.')
    return files_dict



# ----------------- #

def create_matrix_header(files_dict):
    # Merge datasets
    files_merged = reduce(lambda left, right: pd.merge(left, right, on=['phosphosite_ID'], how='outer'), files_dict.values())

    print('Datasets have been merged on phosphosite_ID column.')
    
    # Make sure phosphosites are rows
    matrix_cols = files_merged.set_index('phosphosite_ID').T  # Transpose so phosphosites are rows

    matrix_cols.to_csv('/data/home/bty449/ExplainableAI/RawMatrixProcessing/raw-matrix-header.csv')
    
    print('Unique phosphosite_IDs saved as rows.')
    
    return matrix_cols



# ----------------- #

def add_rows_to_matrix(matrix, files_datasets, files_dict):
    new_columns = []

    for dataset_key, dataset_names in files_datasets:
        if dataset_key in matrix.columns:
            print(f"Dataset key {dataset_key} is already in the matrix.")
            continue

        # Fetch the dataset and set phosphosite_ID as the index
        dataset = files_dict[dataset_key].set_index('phosphosite_ID')
        
        # Make sure each dataset is transposed so that each phosphosite is a row, not a column
        dataset = dataset.T  # Transpose dataset so phosphosites are rows, not columns

        # Rename the index to the dataset name
        dataset.index = [dataset_key]
        
        # Add the dataset as a new column to the matrix (each dataset as a new column)
        matrix = pd.concat([matrix, dataset], axis=0)  # Concatenate along the columns (axis=1)

        print(f"{dataset_key} added to matrix")

    matrix.to_csv('/Users/maryamkoddus/Documents/maryam-ko-QMUL-MSc-Project/01_input_data/MatrixCSVs/intermediary-raw-matrix.csv', index=True)
    print('Intermediary raw matrix saved.')
    
    return matrix



# ----------------- #

def find_outliers_IQR(df): # identify outliers from a dataframe using the Interquartile Range (IQR)
    q1 = df.quantile(0.25)
    q3 = df.quantile(0.75)
    IQR = q3 - q1
    outliers = df[((df < (q1 - 1.5 * IQR)) | (df > (q3 + 1.5 * IQR)))]
    return outliers



# ----------------- #

def drop_outliers_IQR(df): # drop outliers from a dataframe using the Interquartile Range (IQR)
    q1 = df.quantile(0.25)
    q3 = df.quantile(0.75)
    IQR = q3 - q1
    not_outliers = df[~((df < (q1 - 1.5 * IQR)) | (df > (q3 + 1.5 * IQR)))]
    return not_outliers



# ----------------- #

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



# ----------------- #

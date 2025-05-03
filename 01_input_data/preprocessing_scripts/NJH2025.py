import sys
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

dataset = 'NJH2025'

print('Loading raw data for', dataset, '...')
data = pd.read_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/01_input_data/raw_data/40279_2025_2217_MOESM3_ESM.csv', header=0)
print('Raw data loaded.')
data

data.columns = data.columns.str.strip()

# filter data to keep only those with localization probability >= 0.85
data = data[data['Localization prob'] >= 0.85] 

# Filtering out semi-colons from 'Amino acid', 'Positions within proteins', and 'Gene names' columns
data = data[~data['Amino acid'].str.contains(';', na=False)]
data = data[~data['Positions within proteins'].str.contains(';', na=False)]
data = data[~data['Gene names'].str.contains(';', na=False)]
data
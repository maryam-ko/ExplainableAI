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
print(files_dict)  # To verify if 'AMK2021' is in the dictionary

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

AMK2021_conditions = ['Ctrl', '45I', '120R', '30R']
AMK2021_replicates = ['p1', 'p3', 'p5']
AMK2021_names = tuple(f'{condition}_{replicate}' for condition in AMK2021_conditions for replicate in AMK2021_replicates)

CF2017_conditions = ['EOC', 'FTE', 'OSE']
CF2017_replicates = [f'{i+1}' for i in range(4)] + [f'{i+13}' for i in range(5)]  # EOC/FTE (1-4), OSE (26-29)
CF2017_names = tuple(f'{condition} {replicate}' for condition in CF2017_conditions for replicate in CF2017_replicates)

EJN2021_conditions = ['Rest', 'Ex']
EJN2021_treatments = ['Basal', 'Ins']
EJN2021_replicates = [str(i) for i in range(1, 6)]  # Replicates from 1 to 5
EJN2021_names = tuple(f'{replicate}_{condition}_{treatment}' for condition in EJN2021_conditions for treatment in EJN2021_treatments for replicate in EJN2021_replicates)

HS2024_conditions = ['B', 'I', 'PD', 'N']
HS2024_timepoints = ['T1', 'T2', 'T3']
HS2024_replicates = ['N1', 'N2']
HS2024_roots = ['R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'R09', 'R10', 'R11']
HS2024_names = tuple(f'{root}_{condition}_{timepoint}' for root in HS2024_roots for condition in HS2024_conditions for timepoint in HS2024_timepoints)

JC2019_intensity_conditions = ['Normalised intensity (A)', 'Normalised intensity (B)', 'Normalised intensity (C)']
JC2019_samples = ['URO0059', 'URO0126', 'URO0158']
JC2019_root = 'TiOx.Testis'
JC2019_names = tuple(f'{JC2019_root}({sample})_{condition}' for sample in JC2019_samples for condition in JC2019_intensity_conditions)

JJ2018_conditions = ['M/L normalized', 'H/L normalized', 'H/M normalized']
JJ2018_timepoints = ['', '___1', '___2', '___3', ' pP1', ' pP1___1', ' pP1___2', ' pP1___3', 
                     ' pP2', ' pP2___1', ' pP2___2', ' pP2___3', ' pP3', ' pP3___1', ' pP3___2', ' pP3___3']
JJ2018_names = tuple(f'Ratio {condition}{timepoint}' for condition in JJ2018_conditions for timepoint in JJ2018_timepoints)

JVO2010_categories = [
    'Ratio M/L Normalized by Protein Noco45', 'Ratio M/L Normalized by Protein Noco3h', 
    'Ratio M/L Normalized by Protein ThymON', 'Ratio H/M Normalized by Protein Noco45', 
    'Ratio H/M Normalized by Protein Noco3h', 'Ratio H/M Normalized by Protein ThymON', 
    'Ratio Variability [%]', 'Log2 Ratio M/L Normalized aphidicolin', 'Log2 Ratio H/M Normalized aphidicolin'
]
JVO2010_conditions = [
    ' - Log2 mitosis', ' - Log2 G1', ' - Log2 G1/S', ' - Log2 early S', ' - Log2 Late S', ' - Log2 G2'
]
JVO2010_names = tuple(
    f'{category}{condition}' for category in JVO2010_categories for condition in JVO2010_conditions
)

JW2015_intensity_conditions = ['Intensity cap1', 'Intensity cap2', 'Intensity cap3', 'Intensity pre1', 'Intensity pre2', 'Intensity pre3']
JW2015_samples = ['Sample1', 'Sample2', 'Sample3']  # Example sample names
JW2015_replicates = ['Rep1', 'Rep2', 'Rep3']  # Example replicate names
JW2015_names = tuple(f'{sample}_{condition}_{replicate}' for sample in JW2015_samples for condition in JW2015_intensity_conditions for replicate in JW2015_replicates)

KS2014_intensity_conditions = [
    'Intensity', 'Intensity C1', 'Intensity C2', 'Intensity C3', 'Intensity C4', 'Intensity C5', 'Intensity C6', 
    'Intensity E15_1', 'Intensity E15_2', 'Intensity E15_3', 'Intensity E15_4', 'Intensity N1', 'Intensity N2', 
    'Intensity N3', 'Intensity N4', 'Intensity pY_C4', 'Intensity pY_C5', 'Intensity pY_C6', 'Intensity pY_E5_1', 
    'Intensity pY_E5_2', 'Intensity pY_E5_3', 'Intensity pY_E5_4', 'Intensity pY_N1', 'Intensity pY_N2', 'Intensity pY_N3', 
    'Intensity pY_N4', 'Intensity pY_PV1', 'Intensity pY_PV2', 'Intensity pY_PV3', 'Intensity pY_PV4'
]

KS2014_samples = ['SampleA', 'SampleB', 'SampleC']  # Example sample names
KS2014_replicates = ['Rep1', 'Rep2', 'Rep3']  # Example replicate names
KS2014_names = tuple(f'{sample}_{condition}_{replicate}' for sample in KS2014_samples for condition in KS2014_intensity_conditions for replicate in KS2014_replicates)

MV2014_ratio_conditions = [
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
MV2014_samples = ['Sample1', 'Sample2', 'Sample3']  # Example sample names
MV2014_replicates = ['Rep1', 'Rep2', 'Rep3']  # Example replicate names
MV2014_names = tuple(f'{sample}_{condition}_{replicate}' for sample in MV2014_samples for condition in MV2014_ratio_conditions for replicate in MV2014_replicates)

NJH2015_conditions = [
    'Subject 1 115/114 (Exercise/Basal) (normalized)', 
    'Subject 2 117/116 (Exercise/Basal) (normalized)', 
    'Subject 3 116/117 (Exercise/Basal) (normalized)', 
    'Subject 4 114/115 (Exercise/Basal) (normalized)'
]

NJH2015_names = tuple(f'{condition}' for condition in NJH2015_conditions)


NJH2024_subjects = ['Sub1', 'Sub2', 'Sub4', 'Sub6', 'Sub7', 'Sub8', 'Sub9', 'Sub10', 'Sub12', 'Sub13']
NJH2024_conditions = ['MICT_Pre', 'MICT_Mid', 'MICT_Post', 'HIIT_Pre', 'HIIT_Mid', 'HIIT_Post']
NJH2024_names = tuple(f'{subject}_{condition}' for subject in NJH2024_subjects for condition in NJH2024_conditions)

PG2020_conditions = ['CT_A', 'CT_B', 'CT_C', 'CT_D', '2h_A', '2h_B', '2h_C', '2h_D', '4h_A', '4h_B', '4h_C', '4h_D', 
                     '6h_A', '6h_B', '6h_C', '6h_D', '8h_A', '8h_B', '8h_C', '8h_D', '10h_A', '10h_B', '10h_C', '10h_D', 
                     'Mock_10h_A', 'Mock_10h_B', 'Mock_10h_C', 'Mock_10h_D']
PG2020_names = tuple(f'{condition}' for condition in PG2020_conditions)

RAM2015_conditions = ['A', 'B', 'C']
RAM2015_names = tuple(f'Log2 {cond}' for cond in RAM2015_conditions)

RB2022_subjects = [f'Subject{i}' for i in range(1, 9)]
RB2022_exercises = ['Endurance', 'Sprint', 'Strength']
RB2022_timepoints = ['Pre', 'Post', 'Recovery']

RB2022_names = tuple(
    f'{subject}_{timepoint}{exercise}'
    for subject in RB2022_subjects
    for exercise in RB2022_exercises
    for timepoint in RB2022_timepoints
)

RM2009_ratios = ['M/L', 'H/L']
RM2009_timepoints = {
    'M/L': 'M:15 min timepoint;L:0 min. timepoint',
    'H/L': 'H:60 min timepoint;L:0 min. timepoint'
}
RM2009_suffixes = ['', ' Normalized by Corresponding Protein Level']

RM2009_names = tuple(
    f'Ratio {ratio}{suffix} ({RM2009_timepoints[ratio]})'
    for ratio in RM2009_ratios
    for suffix in RM2009_suffixes
)

RN2012_names = tuple(
    [f'H0{i}' for i in range(1, 7)] + [f'L1{i}' for i in range(0, 6)]
)

RNJ2017_conditions = [
    ('Tsup', 'Tstim', ['D1', 'D3']),
    ('Tstim', 'Trest', ['D2', 'D3']),
    ('Tsup', 'Trest', ['D1', 'D2', 'D3']),
]

RNJ2017_names = tuple(
    f'Log2 ratio {cond1}:{cond2}_{day}'
    for cond1, cond2, days in RNJ2017_conditions
    for day in days
)

ZCP2016_ratios = ['H/L']
ZCP2016_replicates = [f'Rep {i}' for i in range(1, 4)]

ZCP2016_names = tuple(
    f'Ratio {ratio} {rep}'
    for ratio in ZCP2016_ratios
    for rep in ZCP2016_replicates
)

ZQ2022_conditions = ['C1', 'C2', 'C3', 'N1', 'N2', 'N3']

ZQ2022_names = tuple(
    f'{condition}' for condition in ZQ2022_conditions
)


# ----------------- #
# PAIR DATASET NAMES WITH FILENAMES
# ----------------- #

files_datasets = [
    ('AST2020', AST2020_names),
    ('AMK2021', AMK2021_names),
    ('CF2017', CF2017_names),
    ('EJN2021', EJN2021_names),
    ('HS2024', HS2024_names),
    ('JC2019', JC2019_names),
    ('JJ2018', JJ2018_names),
    ('JVO2010', JVO2010_names),
    ('JW2015', JW2015_names),
    ('KS2014', KS2014_names),
    ('MV2014', MV2014_names),
    ('NJH2015', NJH2015_names),
    ('NJH2024', NJH2024_names),
    ('PG2020', PG2020_names),
    ('RAM2015', RAM2015_names),
    ('RB2022', RB2022_names),
    ('RM2009', RM2009_names),
    ('RN2012', RN2012_names),
    ('RNJ2017', RNJ2017_names),
    ('ZCP2016', ZCP2016_names),
    ('ZQ2022', ZQ2022_names)
]

intermed_matrix = preprocessing.add_rows_to_matrix(matrix, files_datasets, files_dict)
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

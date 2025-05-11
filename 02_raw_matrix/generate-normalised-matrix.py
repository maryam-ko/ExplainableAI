#!/bin/python


# ----------------- #
# LOAD DEPENDENCIES
# ----------------- #

import sys
import os
import pandas as pd
import numpy as np
from sklearn.preprocessing import MinMaxScaler


grandparent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(grandparent_dir)
from funcs import normalising


# ----------------- #
# IMPORT R-NORMALISED MATRIX 
# ----------------- #

matrix = pd.read_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/RawMatrix_NoOutliers.csv', header = 0)
print(matrix.columns)

matrix = matrix.set_index('DatasetName')

dataset_list = normalising.create_dataframe_per_dataset(matrix)
print(f'Dataframe has been created per dataset')



# ----------------- #
# LOAD DATASET NAMES
# ----------------- #

# define names for each dataset
AST2020_conditions = ['Control', 'anti-CD3', 'anti-CD3+PDL2']
AST2020_replicates = ['', '.1', '.2']
AST2020_names = tuple(f'AST2020_{cond}{rep}' for cond in AST2020_conditions for rep in AST2020_replicates)

AMK2021_conditions = ['Ctrl', '45I', '120R', '30R']
AMK2021_replicates = ['p1', 'p3', 'p5']
AMK2021_names = tuple(f'AMK2021_{condition}_{replicate}' for condition in AMK2021_conditions for replicate in AMK2021_replicates)

CF2017_conditions = ['EOC', 'FTE', 'OSE']
CF2017_replicates = [f'{i+1}' for i in range(4)] + [f'{i+13}' for i in range(5)]  
CF2017_names = tuple(f'CF2017_{condition} {replicate}' for condition in CF2017_conditions for replicate in CF2017_replicates)

CKC2017_conditions = ['0min', '15min', '30min', '60min', '90min', '120min']
CKC2017_replicates = ['p1', 'p2', 'p3']
CKC2017_names = tuple(f'CKC2017_{condition}_{replicate}' for condition in CKC2017_conditions for replicate in CKC2017_replicates)

DB2022_conditions = [
    'L1', 'L1.1', 'L1.2',
    'L2', 'L2.1', 'L2.2',
    'R1', 'R1.1', 'R1.2',
    'R2', 'R2.1', 'R2.2',
    'Intensity_L1', 'Intensity_L2', 'Intensity_R1', 'Intensity_R2'
]

DB2022_names = tuple(f'DB2022_{condition}' for condition in DB2022_conditions)

EJN2021_conditions = ['Rest', 'Ex']
EJN2021_treatments = ['Basal', 'Ins']
EJN2021_replicates = [str(i) for i in range(1, 6)]  
EJN2021_names = tuple(f'EJN2021_{replicate}_{condition}_{treatment}' for condition in EJN2021_conditions for treatment in EJN2021_treatments for replicate in EJN2021_replicates)

FS2019_conditions = ['Ctr', 'dbdb']
FS2019_replicates = ['1', '2', '3']
FS2019_names = tuple(f'FS2019_{condition}{rep}' for condition in FS2019_conditions for rep in FS2019_replicates)

GF2025_samples = [f'P{i}' for i in range(1, 11)]

GF2025_conditions = ['Adh', 'Sph']
GF2025_names = tuple(f'GF2025_{sample}_{condition}' 
                     for sample in GF2025_samples 
                     for condition in GF2025_conditions)

GRW2016_conditions = ['H/L', 'M/L']
GRW2016_groups = ['1', '2', '4', '5a']
GRW2016_replicates = ['1', '2', '3']

GRW2016_names = tuple(f'GRW2016_{condition}_{group}___{replicate}' 
                      for condition in GRW2016_conditions 
                      for group in GRW2016_groups 
                      for replicate in GRW2016_replicates)

HH2022_conditions = [str(i) for i in range(10)]  # '0', '1', '2', ..., '9'
HH2022_replicates = ['p1', 'p2', 'p3']

HH2022_names = tuple(f'HH2022_Reporter_intensity_{condition}_{replicate}' 
                     for condition in HH2022_conditions 
                     for replicate in HH2022_replicates)

HS2024_conditions = ['B', 'I', 'PD', 'N']
HS2024_timepoints = ['T1', 'T2', 'T3']
HS2024_replicates = ['N1', 'N2']
HS2024_roots = ['R03', 'R04', 'R05', 'R06', 'R07', 'R08', 'R09', 'R10', 'R11']
HS2024_names = tuple(f'HS2024_{root}_{condition}_{timepoint}' for root in HS2024_roots for condition in HS2024_conditions for timepoint in HS2024_timepoints)

JAW2011_conditions = ['LS1a', 'LS1b', 'LS2', 'HS']
JAW2011_replicates = ['', '_1', '_2', '_3']

JAW2011_names = tuple(
    f'JAW2011_{cond}{rep}' for cond in JAW2011_conditions for rep in JAW2011_replicates
)

JB2023_intensities = ['Intensity', 'Intensity___1', 'Intensity___2', 'Intensity___3']

JB2023_groups = [
    '20220718_C28798_002_S397439_M2_prim_IIM_R3_Group_4',
    '20220718_C28798_003_S397432_M2_prim_IS_R4_Group_2',
    '20220718_C28798_004_S397437_M2_prim_IIM_R1_Group_4',
    '20220718_C28798_005_S397438_M2_prim_IIM_R2_Group_4',
    '20220718_C28798_007_S397436_M2_prim_FCM_R4_Group_3',
    '20220718_C28798_008_S397429_M2_prim_IS_R1_Group_2',
    '20220718_C28798_009_S397426_M2_prim_NT_R2_Group_1',
    '20220718_C28798_010_S397425_M2_prim_NT_R1_Group_1',
    '20220718_C28798_013_S397433_M2_prim_FCM_R1_Group_3',
    '20220718_C28798_014_S397428_M2_prim_NT_R4_Group_1',
    '20220718_C28798_015_S397430_M2_prim_IS_R2_Group_2',
    '20220718_C28798_016_S397431_M2_prim_IS_R3_Group_2',
    '20220718_C28798_018_S397434_M2_prim_FCM_R2_Group_3',
    '20220718_C28798_019_S397440_M2_prim_IIM_R4_Group_4',
    '20220718_C28798_020_S397435_M2_prim_FCM_R3_Group_3',
    '20220718_C28798_021_S397427_M2_prim_NT_R3_Group_1'
]

JB2023_replicates = ['1', '2', '3']

JB2023_names = tuple(f'JB2023_{intensity}_{group}___{replicate}' 
                     for intensity in JB2023_intensities 
                     for group in JB2023_groups 
                     for replicate in JB2023_replicates)

JC2019_intensity_conditions = ['Normalised intensity (A)', 'Normalised intensity (B)', 'Normalised intensity (C)']
JC2019_samples = ['URO0059', 'URO0126', 'URO0158']
JC2019_root = 'TiOx.Testis'
JC2019_names = tuple(f'JC2019_{JC2019_root}({sample})_{condition}' for sample in JC2019_samples for condition in JC2019_intensity_conditions)

JJ2018_conditions = ['M/L normalized', 'H/L normalized', 'H/M normalized']
JJ2018_timepoints = ['', '___1', '___2', '___3', ' pP1', ' pP1___1', ' pP1___2', ' pP1___3', 
                     ' pP2', ' pP2___1', ' pP2___2', ' pP2___3', ' pP3', ' pP3___1', ' pP3___2', ' pP3___3']
JJ2018_names = tuple(f'JJ2018_{condition}{timepoint}' for condition in JJ2018_conditions for timepoint in JJ2018_timepoints)

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
    f'JVO2010_{category}{condition}' for category in JVO2010_categories for condition in JVO2010_conditions
)

JW2015_intensity_conditions = ['Intensity cap1', 'Intensity cap2', 'Intensity cap3', 'Intensity pre1', 'Intensity pre2', 'Intensity pre3']
JW2015_samples = ['Sample1', 'Sample2', 'Sample3'] 
JW2015_replicates = ['Rep1', 'Rep2', 'Rep3']  
JW2015_names = tuple(f'JW2015_{sample}_{condition}_{replicate}' for sample in JW2015_samples for condition in JW2015_intensity_conditions for replicate in JW2015_replicates)

LATM2013_conditions = ['Forward', 'Reverse']
LATM2013_types = ['Raw', 'Normalized']

LATM2013_names = tuple(
    f'LATM2013_{t}_{c}' for t in LATM2013_types for c in LATM2013_conditions
)

LG2023_conditions = ['GSK3α', 'Control']
LG2023_intensities = ['1', '2', '3']
LG2023_replicates = ['1', '2', '3']

LG2023_names = tuple(f'LG2023_{condition}-{intensity}_intensity_{replicate}' 
                     for condition in LG2023_conditions 
                     for intensity in LG2023_intensities 
                     for replicate in LG2023_replicates)

KS2014_intensity_conditions = [
    'Intensity', 'Intensity C1', 'Intensity C2', 'Intensity C3', 'Intensity C4', 'Intensity C5', 'Intensity C6', 
    'Intensity E15_1', 'Intensity E15_2', 'Intensity E15_3', 'Intensity E15_4', 'Intensity N1', 'Intensity N2', 
    'Intensity N3', 'Intensity N4', 'Intensity pY_C4', 'Intensity pY_C5', 'Intensity pY_C6', 'Intensity pY_E5_1', 
    'Intensity pY_E5_2', 'Intensity pY_E5_3', 'Intensity pY_E5_4', 'Intensity pY_N1', 'Intensity pY_N2', 'Intensity pY_N3', 
    'Intensity pY_N4', 'Intensity pY_PV1', 'Intensity pY_PV2', 'Intensity pY_PV3', 'Intensity pY_PV4'
]

KS2014_samples = ['SampleA', 'SampleB', 'SampleC']  
KS2014_replicates = ['Rep1', 'Rep2', 'Rep3'] 
KS2014_names = tuple(f'KS2014_{sample}_{condition}_{replicate}' for sample in KS2014_samples for condition in KS2014_intensity_conditions for replicate in KS2014_replicates)

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
MV2014_samples = ['Sample1', 'Sample2', 'Sample3']  
MV2014_replicates = ['Rep1', 'Rep2', 'Rep3'] 
MV2014_names = tuple(f'MV2014_{sample}_{condition}_{replicate}' for sample in MV2014_samples for condition in MV2014_ratio_conditions for replicate in MV2014_replicates)

NJH2015_conditions = [
    'Subject 1 115/114 (Exercise/Basal) (normalized)', 
    'Subject 2 117/116 (Exercise/Basal) (normalized)', 
    'Subject 3 116/117 (Exercise/Basal) (normalized)', 
    'Subject 4 114/115 (Exercise/Basal) (normalized)'
]

NJH2015_names = tuple(f'NJH2015_{condition}' for condition in NJH2015_conditions)


NJH2024_subjects = ['Sub1', 'Sub2', 'Sub4', 'Sub6', 'Sub7', 'Sub8', 'Sub9', 'Sub10', 'Sub12', 'Sub13']
NJH2024_conditions = ['MICT_Pre', 'MICT_Mid', 'MICT_Post', 'HIIT_Pre', 'HIIT_Mid', 'HIIT_Post']
NJH2024_names = tuple(f'NJH2024_{subject}_{condition}' for subject in NJH2024_subjects for condition in NJH2024_conditions)

NJH2025_subjects = ['Sub1', 'Sub2', 'Sub4', 'Sub6', 'Sub7', 'Sub8', 'Sub9', 'Sub10', 'Sub12', 'Sub13']

NJH2025_exercise_types = ['MICT', 'HIIT']
NJH2025_time_points = ['Pre', 'Mid', 'Post']
NJH2025_names = tuple(f'NJH2025_{subject}_{exercise_type}_{time_point}' 
                      for subject in NJH2025_subjects
                      for exercise_type in NJH2025_exercise_types
                      for time_point in NJH2025_time_points)

PG2020_conditions = ['CT_A', 'CT_B', 'CT_C', 'CT_D', '2h_A', '2h_B', '2h_C', '2h_D', '4h_A', '4h_B', '4h_C', '4h_D', 
                     '6h_A', '6h_B', '6h_C', '6h_D', '8h_A', '8h_B', '8h_C', '8h_D', '10h_A', '10h_B', '10h_C', '10h_D', 
                     'Mock_10h_A', 'Mock_10h_B', 'Mock_10h_C', 'Mock_10h_D']
PG2020_names = tuple(f'PG2020_{condition}' for condition in PG2020_conditions)

RAM2015_conditions = ['A', 'B', 'C']
RAM2015_names = tuple(f'RAM2015_Log2 {cond}' for cond in RAM2015_conditions)

RB2022_subjects = [f'Subject{i}' for i in range(1, 9)]
RB2022_exercises = ['Endurance', 'Sprint', 'Strength']
RB2022_timepoints = ['Pre', 'Post', 'Recovery']

RB2022_names = tuple(
    f'RB2022_{subject}_{timepoint}{exercise}'
    for subject in RB2022_subjects
    for exercise in RB2022_exercises
    for timepoint in RB2022_timepoints
)

RBK2021_time_points = ['1244T', '1259T', '1260T', '1326T', '1332T']

RBK2021_phases = ['Early', 'Late']

RBK2021_names = tuple(f'RBK2021_{time_point}_{phase}' 
                      for time_point in RBK2021_time_points
                      for phase in RBK2021_phases)


RM2009_ratios = ['M/L', 'H/L']
RM2009_timepoints = {
    'M/L': 'M:15 min timepoint;L:0 min. timepoint',
    'H/L': 'H:60 min timepoint;L:0 min. timepoint'
}
RM2009_suffixes = ['', ' Normalized by Corresponding Protein Level']

RM2009_names = tuple(
    f'RM2009_Ratio {ratio}{suffix} ({RM2009_timepoints[ratio]})'
    for ratio in RM2009_ratios
    for suffix in RM2009_suffixes
)

RN2012_names = tuple(
    [f'RN2012_H0{i}' for i in range(1, 7)] + [f'RN2012_L1{i}' for i in range(0, 6)]
)

RNJ2017_conditions = [
    ('Tsup', 'Tstim', ['D1', 'D3']),
    ('Tstim', 'Trest', ['D2', 'D3']),
    ('Tsup', 'Trest', ['D1', 'D2', 'D3']),
]

RNJ2017_names = tuple(
    f'RNJ2017_Log2 ratio {cond1}:{cond2}_{day}'
    for cond1, cond2, days in RNJ2017_conditions
    for day in days
)

SFR2015_samples = ['SampleA', 'SampleB', 'SampleC']
SFR2015_conditions = ['a', 'b', 'c']  
SFR2015_replicates = ['Rep1', 'Rep2', 'Rep3']

SFR2015_names = tuple(f'SFR2015_{sample}_HL_log2_ratio_avg_{condition}_{replicate}' 
                      for sample in SFR2015_samples 
                      for condition in SFR2015_conditions 
                      for replicate in SFR2015_replicates)

TM2022_treatment_types = ['S2.4h', 'S2.8h', 'S2.12h', 'mock']

TM2022_time_points = ['A', 'B', 'C']

TM2022_names = tuple(f'TM2022_{treatment}_{time_point}' 
                     for treatment in TM2022_treatment_types
                     for time_point in TM2022_time_points)

YB2013_experiments = ['experiment 1', 'experiment 2']

YB2013_measurement = 'Ratio H/L Normalized'

YB2013_names = tuple(f'YB2013_{measurement}_in_{experiment}' 
                     for experiment in YB2013_experiments
                     for measurement in [YB2013_measurement])

ZCP2016_ratios = ['H/L']
ZCP2016_replicates = [f'ZCP2016_Rep {i}' for i in range(1, 4)]

ZCP2016_names = tuple(
    f'ZCP2016_Ratio {ratio} {rep}'
    for ratio in ZCP2016_ratios
    for rep in ZCP2016_replicates
)

ZQ2022_conditions = ['C1', 'C2', 'C3', 'N1', 'N2', 'N3']

ZQ2022_names = tuple(
    f'ZQ2022_{condition}' for condition in ZQ2022_conditions
)

all_names = (AST2020_names + AMK2021_names + CF2017_names + CKC2017_names + DB2022_names + EJN2021_names +
            FS2019_names + GF2025_names + GRW2016_names + HH2022_names + HS2024_names + JAW2011_names +
            JB2023_names + JC2019_names + JJ2018_names + JVO2010_names + JW2015_names + LATM2013_names +
            LG2023_names + KS2014_names + MV2014_names + NJH2015_names + NJH2024_names + NJH2025_names +
            PG2020_names + RAM2015_names + RB2022_names + RBK2021_names + RM2009_names + RN2012_names +
            RNJ2017_names + SFR2015_names + TM2022_names + YB2013_names + ZCP2016_names + ZQ2022_names)


# ----------------- #
# CREATE DICTIONARY OF DATAFRAMES
# ----------------- #

dataset_dict = {}

# add each dataframe from list to dictionary (with index)
for i, dataset_name in enumerate(all_names):
    dataset_dict[dataset_name] = dataset_list[i]

print(dataset_dict[f'AST2020_Control'].head())  


# ----------------- #
# CREATE MATRIX HEADER
# ----------------- #

norm_matrix = normalising.MinMax_normalize_and_merge(dataset_dict, MinMaxScaler())

norm_matrix.columns.name = ''

# get unique phosphosites
uniq_phos = norm_matrix['DatasetName'].unique()

if len(uniq_phos) == len(norm_matrix):
    print('Passed! Number of unique phosphosites matches length of merged dataframe')
else:
    print('Failed! Number of unique phosphosites does not match length of merged dataframe')

matrix_cols = pd.DataFrame(columns = uniq_phos)

matrix_cols.to_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/normalised-matrix-header.csv', index = False)

print('Matrix header saved to CSV.', matrix_cols)


# ----------------- #
# FORMAT MATRIX
# ----------------- #

transposed_matrix = norm_matrix.T
print(f'Transposed matrix:', transposed_matrix.head())

# set first row to column header
transposed_matrix.columns = transposed_matrix.iloc[0]

# remove first row
transposed_matrix = transposed_matrix[1:]

# rename and reset column index
transposed_matrix = transposed_matrix.rename_axis('DatasetName').reset_index()
print(f'Transposed matrix:', transposed_matrix.head())


cols = [i for i in transposed_matrix.columns if i not in ['DatasetName']]

for col in cols:
    transposed_matrix[col] = pd.to_numeric(transposed_matrix[col])
    
minmax_matrix = transposed_matrix.replace([np.inf, -np.inf], np.nan)
print(f'MinMax matrix:', minmax_matrix.head())

minmax_matrix.to_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/NormalisedMatrix.csv', index = False)
print(f'Normalised matrix saved successfully!', minmax_matrix)
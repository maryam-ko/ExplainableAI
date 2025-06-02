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
matrix = pd.read_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NBA-Matrix_Quantile.csv', header = 0)

matrix = matrix.set_index('DatasetName')

dataset_list = normalising.create_dataframe_per_dataset(matrix)
print(f'Dataframe has been created per dataset')


# ----------------- #
# LOAD DATASET NAMES
# ----------------- #

# define names for each dataset
# AB2020_names = (
#     [f'AB2020_Intensity___{i}' for i in range(1, 4)] +
#     [f'AB2020_Intensity {cond}_{rep}' for cond in ['EMD', 'EMD_IR', 'IR', 'untreated'] for rep in range(1, 3)] +
#     [f'AB2020_Intensity {cond}_{rep}___{r}' for cond in ['EMD', 'EMD_IR', 'IR', 'untreated'] for rep in range(1, 3) for r in range(1, 4)]
# )

AH2018_names = (
    [f'AH2018_TMTMS{i}_{rep}_{cond}' for i in [2, 3] for rep in range(1, 4) for cond in ['C', 'X', 'N']] +
    [f'AH2018_TMTMS{i}_Pool' for i in [2, 3]] +
    [f'AH2018_LFQ_{rep}_{cond}' for rep in range(1, 5) for cond in ['C', 'N', 'X']] +
    [f'AH2018_LFQ_M_{rep}_{cond}' for rep in range(1, 5) for cond in ['C', 'N', 'X']] +
    [f'AH2018_LFQ_IT_{rep}_{cond}' for rep in [1, 3, 4] for cond in ['C', 'X', 'N']] +
    [f'AH2018_SILAC_{rep}_{cond}' for rep in range(1, 4) for cond in ['C', 'X', 'N']] +
    [f'AH2018_SILAC_M_{rep}_{cond}' for rep in range(1, 4) for cond in ['C', 'X', 'N']] +
    [f'AH2018_SILAC_MR_{rep}_{cond}' for rep in range(1, 4) for cond in ['C', 'X', 'N']] +
    [f'AH2018_SILAC_IT_{rep}_{cond}' for rep in range(1, 4) for cond in ['C', 'X', 'N']]
)


AMK2021_names = (
    [f'AMK2021_{cond}_p1' for cond in ['Ctrl', '45I', '120R']] +
    [f'AMK2021_{cond}_p3' for cond in ['Ctrl', '45I', '30R', '120R']] +
    [f'AMK2021_{cond}_p5' for cond in ['Ctrl', '45I', '30R', '120R']]
)


AR2021_names = (
    [f'AR2021_L41 {dose}' for dose in ['1 μM', '10 μM']] +
    [f'AR2021_ALG {dose}' for dose in ['1 μM', '10 μM']] +
    [f'AR2021_KD {time}' for time in ['48h', '72h']]
)


AS2011_names = [
    'AS2011_Inhi 1_Ratio H/L normalized by corresponding protein level (H:TAL treatment; L: MA treatment)',
    'AS2011_Inhi 2_Ratio H/L normalized by corresponding protein level (H:TAL treatment; L: MA treatment)',
    'AS2011_RNAi 1_Ratio H/L normalized by corresponding protein level (H: shPlk1; L: shEg5)',
    'AS2011_RNAi 2_Ratio H/L normalized by corresponding protein level (H: shPlk1; L: shEg5)',
    
    'AS2011_Ratio H/L',
    *[f'AS2011_Ratio H/L_{i}' for i in range(1, 4)],
    
    'AS2011_Ratio H/L Normalized',
    *[f'AS2011_Ratio H/L Normalized_{i}' for i in range(1, 4)],
    
    'AS2011_Ratio H/L Normalized by Proteins',
    'AS2011_Ratio H/L Normalized by Protein',
    'AS2011_Ratio H/L Normalized by Protein (Unnormalized)',
    
    'AS2011_Ratio H/L Proteins',
    'AS2011_Ratio H/L Protein',
]

conditions = ['Inhi_1', 'Inhi_2', 'RNAi_1', 'RNAi_2']
prefixes = [
    'Ratio H/L', 'Ratio H/L Normalized', 
    'Ratio H/L Normalized by Proteins', 'Ratio H/L Normalized by Protein',
    'Ratio H/L Proteins', 'Ratio H/L Protein'
]

for cond in conditions:
    AS2011_names += [f'AS2011_{prefix} {cond}' for prefix in prefixes]
    AS2011_names += [f'AS2011_Ratio H/L {cond}_{i}' for i in range(1, 4)]
    AS2011_names += [f'AS2011_Ratio H/L Normalized {cond}_{i}' for i in range(1, 4)]



AST2020_names = (
    [f'AST2020_{condition}' for condition in ['Control', 'anti-CD3', 'anti-CD3+PDL2']] +
    [f'AST2020_{condition}{replicate}' for condition in ['Control', 'anti-CD3', 'anti-CD3+PDL2'] for replicate in ['.1', '.2']]
)


BS2014_names = (
    [f'BS2014_Log2(Ratio H/L normalized Rep_{i})' for i in range(1, 4)] + 
    ['BS2014_Log2(Ratio H/L normalized)', 'BS2014_Log2(Ratio H/L unmodified peptide)', 'BS2014_Log2(Ratio H/L protein)']
)


CF2017_names = (
    [f'CF2017_EOC {i}' for i in range(1, 5)] +
    [f'CF2017_FTE {i}' for i in range(13, 18)] +
    [f'CF2017_OSE {i}' for i in range(26, 30)]
)
                 
CK2022_names = (
    [f'CK2022_{i} 16-1' for i in range(1, 17)] +
    [f'CK2022_{i} 200-1' for i in range(1, 17)] +
    [f'CK2022_{i} 40-1' for i in range(1, 17)]
)

CKC2017_names = (
    [f'CKC2017_log2 normalized phosphoproteomics_{time}min_rep#{rep}' for time in [0, 15, 30, 60, 90, 120] for rep in range(1, 4)]
)

CP2009_names = [
    'CP2009_Ratio M/L',
    'CP2009_Ratio M/L Normalized (5 nM dasatinib vs Ctrl)',
    'CP2009_Ratio H/L',
    'CP2009_Ratio H/L Normalized (50 nM dasatinib vs Ctrl)',
    'CP2009_Ratio H/M',
    'CP2009_Ratio H/M Normalized (50 nM dasatinib vs 5 nM dasatinib)'
]

# CR2017_names = (
#     [f'CR2017_Ratio Replicate {i}' for i in range(1, 4)] +
#     [f'CR2017_Ratio Replicate {i}.{j}' for j in range(1, 6) for i in range(1, 4)]
# )

CW2012_names = (
    [f'CW2012_10 µM Erl/control Experiment {exp} Replicate {rep}' for exp in range(1, 4) for rep in ['A', 'B']] +
    [f'CW2012_10 µM Gef/control Experiment {exp} Replicate {rep}' for exp in range(1, 4) for rep in ['A', 'B']]
)

DB2022_names = (
    [f'DB2022_LFQ intensity {side}{rep}' for side in ['L1', 'L2', 'R1', 'R2'] for rep in ['', '.1', '.2']] +
    [f'DB2022_Intensity {side}' for side in ['L1', 'L2', 'R1', 'R2']]
)

DL2023_names = (
    [f'DL2023_control-{i}' for i in range(1, 4)] +
    [f'DL2023_12h-{i}' for i in range(1, 4)] +
    [f'DL2023_24h-{i}' for i in range(1, 4)]
)

EJN2021_names = (
    [f'EJN2021_{i}_Rest_Basal' for i in range(1, 6)] +
    [f'EJN2021_{i}_Ex_Basal' for i in range(1, 6)] +
    [f'EJN2021_{i}_Rest_Ins' for i in range(1, 6)] +
    [f'EJN2021_{i}_Ex_Ins' for i in range(1, 6)]
)

# EP2017_names = (
#     [f'EP2017_Untreated R#{i} log2 ratio' for i in range(1, 5)] +
#     [f'EP2017_Pervanadate R#{i} log2 ratio' for i in range(1, 4)] +
#     [f'EP2017_Mitotic R#{i} log2 ratio' for i in range(1, 4)] +
#     [f'EP2017_Untreated R#{i} normalized log2 ratio' for i in range(1, 5)] +
#     [f'EP2017_Pervanadate R#{i} normalized log2 ratio' for i in range(1, 4)] +
#     [f'EP2017_Mitotic R#{i} normalized log2 ratio' for i in range(1, 4)]
# )

FR2021_names = (
    [f'FR2021_Intensity.{label}' for label in ['A', 'B', 'G', 'H', 'M', 'N']] +
    [f'FR2021_Norm.Intensity.{label}' for label in ['A', 'B', 'G', 'H', 'M', 'N']]
)

FS2019_names = (
    [f'FS2019_Ctr{i}' for i in range(1, 4)] +
    [f'FS2019_db/db{i}' for i in range(1, 4)]
)

FSO2012_names = (
    [f'FSO2012_Ratio -/+ 3-MB-PP1 exp. as{i}{suffix}' for i in range(1, 5) for suffix in ['a', 'b']] +
    [f'FSO2012_Ratio -/+ 3-MB-PP1 exp. wt{i}{suffix}' for i in range(1, 5) for suffix in ['a', 'b']] +
    [f'FSO2012_Ratio -/+ 3-MB-PP1 exp. as{i}' for i in range(1, 5)] +
    [f'FSO2012_Ratio -/+ 3-MB-PP1 exp. wt{i}' for i in range(1, 5)]
)

FSO2013_names = (
    [f'FSO2013_Ratio M/L normalized Rep{i}{suffix}' for i in range(1, 4) for suffix in ['', '___1', '___2', '___3']] +
    [f'FSO2013_Ratio H/L normalized Rep{i}{suffix}' for i in range(1, 4) for suffix in ['', '___1', '___2', '___3']] +
    [f'FSO2013_Ratio H/M normalized Rep{i}{suffix}' for i in range(1, 4) for suffix in ['', '___1', '___2', '___3']]
)

GF2025_names = (
    [f'GF2025_P{i}_{condition}' for i in range(1, 11) for condition in ['Adh', 'Sph']]
)

GRW2016_names = [
    f'GRW2016_Ratio H/L normalized SILAC {i}___{j}' for i in [1, 2, 4] for j in range(1, 4)
] + [
    f'GRW2016_Ratio M/L normalized SILAC 5a___{i}' for i in range(1, 4)
]

HH2022_names = [
    f'HH2022_Reporter intensity {i}' for i in range(10)
]

HS2024_names = (
    [f'HS2024_{i}_{condition}' for i in range(3, 12) for condition in ['B_T1', 'B_T2', 'B_T3', 'I_T1', 'I_T2', 'I_T3', 'PD_T1', 'PD_T2', 'PD_T3', 'B_N1', 'B_N2', 'I_N1', 'I_N2']]
)

HW2015_names = (
    [f'HW2015_Ratio M/L normalized Exp{i}' for i in range(1, 4)] +
    [f'HW2015_Ratio H/L normalized Exp{i}' for i in range(1, 4)]
)

IMM2018_names = (
    [f'IMM2018_normalised rep. int. CDK5L-{i}' for i in range(1, 4)] +
    [f'IMM2018_normalised rep. int. Empty-{i}' for i in range(1, 4)]
)

# IUA2019_names = (
#     [f'IUA2019_CONTROL_{i} Soluble_1min' for i in range(1, 4)] +
#     [f'IUA2019_U50488H_{i} Soluble_1min' for i in range(1, 4)] +
#     [f'IUA2019_CONTROL_{i} Soluble_60minutes' for i in range(1, 4)] +
#     [f'IUA2019_U50488H_{i} Soluble_60minutes' for i in range(1, 4)] +
#     [f'IUA2019_CONTROL_{i} Insoluble_1min' for i in range(1, 4)] +
#     [f'IUA2019_U50488H_{i} Insoluble_1min' for i in range(1, 4)] +
#     [f'IUA2019_CONTROL_{i} Insoluble_60minutes' for i in range(1, 4)] +
#     [f'IUA2019_U50488H_{i} Insoluble_60minutes' for i in range(1, 4)]
# )

JAW2011_names = (
    [f'JAW2011_Ratio H/L'] +
    [f'JAW2011_Ratio H/L_{i}' for i in range(1, 4)] +
    [f'JAW2011_Ratio H/L LS1a'] +
    [f'JAW2011_Ratio H/L LS1a_{i}' for i in range(1, 4)] +
    [f'JAW2011_Ratio H/L LS1b'] +
    [f'JAW2011_Ratio H/L LS1b_{i}' for i in range(1, 4)] +
    [f'JAW2011_Ratio H/L LS2'] +
    [f'JAW2011_Ratio H/L LS2_{i}' for i in range(1, 4)] +
    [f'JAW2011_Ratio H/L HS'] +
    [f'JAW2011_Ratio H/L HS_{i}' for i in range(1, 4)]
)

# groups = [
#     '20220718_C28798_002_S397439_M2_prim_IIM_R3_Group_4', '20220718_C28798_003_S397432_M2_prim_IS_R4_Group_2',
#     '20220718_C28798_004_S397437_M2_prim_IIM_R1_Group_4', '20220718_C28798_005_S397438_M2_prim_IIM_R2_Group_4',
#     '20220718_C28798_007_S397436_M2_prim_FCM_R4_Group_3', '20220718_C28798_008_S397429_M2_prim_IS_R1_Group_2',
#     '20220718_C28798_009_S397426_M2_prim_NT_R2_Group_1', '20220718_C28798_010_S397425_M2_prim_NT_R1_Group_1',
#     '20220718_C28798_013_S397433_M2_prim_FCM_R1_Group_3', '20220718_C28798_014_S397428_M2_prim_NT_R4_Group_1',
#     '20220718_C28798_015_S397430_M2_prim_IS_R2_Group_2', '20220718_C28798_016_S397431_M2_prim_IS_R3_Group_2',
#     '20220718_C28798_018_S397434_M2_prim_FCM_R2_Group_3', '20220718_C28798_019_S397440_M2_prim_IIM_R4_Group_4',
#     '20220718_C28798_020_S397435_M2_prim_FCM_R3_Group_3', '20220718_C28798_021_S397427_M2_prim_NT_R3_Group_1'
# ]

# JB2023_names = (
#     [f'JB2023_Intensity'] +
#     [f'JB2023_Intensity___{i}' for i in range(1, 4)] +
#     [f'JB2023_Intensity {group}' for group in groups] +
#     [f'JB2023_Intensity {group}___{i}' for group in groups for i in range(1, 4)]
# )

JC2019_names = [
    'JC2019_Normalised_intensity_A_TiOx_Testis_sample_URO0059',
    'JC2019_Normalised_intensity_B_TiOx_Testis_sample_URO0126',
    'JC2019_Normalised_intensity_C_TiOx_Testis_sample_URO0158'
]

JJ2018_names = (
    [f'JJ2018_Ratio M/L normalized{suffix}' for suffix in ['', '___1', '___2', '___3']] +
    [f'JJ2018_Ratio M/L normalized pP{i}{suffix}' for i in range(1, 4) for suffix in ['', '___1', '___2', '___3']] +
    [f'JJ2018_Ratio H/L normalized{suffix}' for suffix in ['', '___1', '___2', '___3']] +
    [f'JJ2018_Ratio H/L normalized pP{i}{suffix}' for i in range(1, 4) for suffix in ['', '___1', '___2', '___3']] +
    [f'JJ2018_Ratio H/M normalized{suffix}' for suffix in ['', '___1', '___2', '___3']] +
    [f'JJ2018_Ratio H/M normalized pP{i}{suffix}' for i in range(1, 4) for suffix in ['', '___1', '___2', '___3']]
)

JVO2010_names = (
    [f'JVO2010_Ratio M/L Normalized by Protein Noco45 - Log2 mitosis'] +
    [f'JVO2010_Ratio M/L Normalized by Protein Noco3h - Log2 G1'] +
    [f'JVO2010_Ratio M/L Normalized by Protein ThymON - Log2 G1/S'] +
    [f'JVO2010_Ratio H/M Normalized by Protein Noco45 - Log2 early S'] +
    [f'JVO2010_Ratio H/M Normalized by Protein Noco3h - Log2 Late S'] +
    [f'JVO2010_Ratio H/M Normalized by Protein ThymON - Log2 G2'] +
    [f'JVO2010_Log2 Ratio M/L Normalized aphidicolin'] +
    [f'JVO2010_Log2 Ratio H/M Normalized aphidicolin']
)

JW2015_names = (
    [f'JW2015_Intensity cap{i}' for i in range(1, 4)] +
    [f'JW2015_Intensity pre{i}' for i in range(1, 4)]
)

JW2021_names = (
    [f'JW2021_Before_{i}' for i in range(1, 4)] +  
    [f'JW2021_After_{i}' for i in range(1, 4)]    
)

# KBE2022_names = [
#     f'KBE2022_Log2_Ex_{i}_{suffix}' 
#     for i in [13, 14, 15, 18, 19, 24, 25, 27, 28, 29, 34, 36, 39, 40, 41, 42, 44, 6, 8, 9] 
#     for suffix in ['T', 'U']
# ]

# KDH2022_names = [
#     f'KDH2022_Intensity_{label}' 
#     for label in ['D1', 'D2', 'D3', 'D4', 'D5', 'N1', 'N2', 'N3', 'N4', 'N5']
# ] 

KS2014_names = [
    f'KS2014_Intensity{suffix}' for suffix in [
        '', ' C1', ' C2', ' C3', ' C4', ' C5', ' C6', 
        ' E15_1', ' E15_2', ' E15_3', ' E15_4', 
        ' N1', ' N2', ' N3', ' N4', 
        ' pY_C4', ' pY_C5', ' pY_C6', 
        ' pY_E5_1', ' pY_E5_2', ' pY_E5_3', ' pY_E5_4', 
        ' pY_N1', ' pY_N2', ' pY_N3', ' pY_N4', 
        ' pY_PV1', ' pY_PV2', ' pY_PV3', ' pY_PV4'
    ]
]

# KTGR2014_names = (
#     [f'KTGR2014_Ratio-AA {time} mins' for time in [2, 7, 15, 30]] +
#     [f'KTGR2014_Ratio-Rapa {time} mins' for time in [2, 7, 15, 30]] +
#     [f'KTGR2014_Ratio M L Normalized Rep{i}_{t}' for i in range(1, 3) for t in ['AA_starvation_0_15', 'Rapamycin_0_15']] +
#     [f'KTGR2014_Ratio H L Normalized Rep{i}_{t}' for i in range(1, 3) for t in ['AA_starvation_0_15', 'Rapamycin_0_15']] +
#     [f'KTGR2014_Ratio H M Normalized Rep{i}_{t}' for i in range(1, 3) for t in ['AA_starvation_0_15', 'Rapamycin_0_15']]
# )

LATM2013_names = [
    'LATM2013_Log2 Ratio (Medium/Light) Forward',
    'LATM2013_Log2 Ratio (Light/Medium) Reverse',
    'LATM2013_Normalized Log2 Ratio (Medium/Light) Forward',
    'LATM2013_Normalized Log2 Ratio (Light/Medium) Reverse',
]

LG2023_names = [
    f'LG2023_GSK3α-{i}_intensity' for i in range(1, 4)
] + [
    f'LG2023_Control-{i}_intensity' for i in range(1, 4)
]

LY2019_names = [
    'LY2019_P1+', 'LY2019_P2+', 'LY2019_P3+', 'LY2019_P4+',
    'LY2019_P1-', 'LY2019_P2-', 'LY2019_P3-', 'LY2019_P4-',
    'LY2019_P5-', 'LY2019_P6-', 'LY2019_P7-', 'LY2019_P8-',
    'LY2019_P5+', 'LY2019_P6+', 'LY2019_P7+', 'LY2019_P8+'
]

MM2018_names = [
    'MM2018_Ratio H/L replicate1',
    'MM2018_Ratio H/L replicate2',
    'MM2018_Ratio H/L replicate 3'
]

MV2014_names = (
    [f'MV2014_Ratio M/L normalized A1___{i} (2 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio H/M normalized A2___{i} (2 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio L/H normalized A3___{i} (2 min)' for i in range(1, 4)] +

    [f'MV2014_Ratio H/L normalized A1___{i} (10 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio L/M normalized A2___{i} (10 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio M/H normalized A3___{i} (10 min)' for i in range(1, 4)] +

    [f'MV2014_Ratio M/L normalized B1___{i} (5 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio H/M normalized B2___{i} (5 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio L/H normalized B3___{i} (5 min)' for i in range(1, 4)] +

    [f'MV2014_Ratio H/L normalized B1___{i} (30 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio L/M normalized B2___{i} (30 min)' for i in range(1, 4)] +
    [f'MV2014_Ratio M/H normalized B3___{i} (30 min)' for i in range(1, 4)]
)

NJH2015_names = (
    [f'NJH2015_iTRAQ Subject {i} {label} (Exercise/Basal) (normalized)' for i, label in zip(range(1, 5), ['115/114', '117/116', '116/117', '114/115'])] +
    [f'NJH2015_TMT Subject {i} {label} (Exercise/Basal) (normalized)' for i, label in zip(range(1, 5), ['127_N/126', '128_N/127_C', '129_N/128_C', '130_N/129_C'])] +
    ['NJH2015_TMT Pool 131/130_C (Exercise/Basal) (normalized)']
)

NJH2024_names = (
    [f'NJH2024_Sub{i}_{exercise}_{time}' for i in [1, 2, 4, 6, 7, 8, 9, 10, 12, 13] for exercise in ['MICT', 'HIIT'] for time in ['Pre', 'Mid', 'Post']]
)

NJH2025_names = (
    [f'NJH2025_Sub{i}_{exercise}_{time}' for i in [1, 2, 4, 6, 7, 8, 9, 10, 12, 13] for exercise in ['MICT', 'HIIT'] for time in ['Pre', 'Mid', 'Post']]
)

PG2020_names = (
    [f'PG2020_Intensity {time}_{condition}' for time in ['CT', '2h', '4h', '6h', '8h', '10h', 'Mock_10h']
     for condition in ['A', 'B', 'C', 'D']]
)

RAM2015_names = ['RAM2015_Log2 A', 'RAM2015_Log2 B', 'RAM2015_Log2 C']

RB2022_names = [
    f'RB2022_{i}_{condition}'
    for i in range(1, 9)  
    for condition in ['PreEndurance', 'PostEndurance', 'RecoveryEndurance', 
                      'PreSprint', 'PostSprint', 'RecoverySprint', 
                      'PreStrength', 'PostStrength', 'RecoveryStrength']
]

RBK2021_names = [
    f'RBK2021_{subject} {phase}' 
    for subject, phase in zip(['1244T', '1259T', '1260T', '1326T', '1332T'], ['Early', 'Early', 'Early', 'Late', 'Late'])
]

RCJD2014_names = [
    f'RCJD2014_{time}_1' for time in ['5min', '10min', '15min', '20min']
] + [
    f'RCJD2014_{time}_2' for time in ['5min', '10min', '15min', '20min']
] + [
    f'RCJD2014_{time}_3' for time in ['5min', '10min', '15min', '20min']
] + [
    f'RCJD2014_{time}_mean' for time in ['5min', '10min', '15min', '20min']
]

RKK2015_names = [
    f'RKK2015_{condition}' for condition in [
        'TNF-α/basal', 'SC-514/basal', 'TNF-α/SC-514', 'TNF-α/basal.1', 'TNF-α+SC-514/basal', 
        'TNF-α/TNF-α+SC-514', 'TNF-α/IKKb', 'IKK+TNF-α/IKKb', 'TNF-α+IKKb/TNF-α', 
        'TNF-α/IKKb(K44M)', 'IKKb(K44M)+TNF-α/IKKb(K44M)', 'IKKb(K44M)+TNF-α/TNF-α'
    ]
]

RM2009_names = [
    f'RM2009_{name}' for name in [
        'Ratio M/L (M:15 min timepoint;L:0 min. timepoint)',
        'Ratio M/L Normalized by Corresponding Protein Level (M:15 min timepoint;L:0 min. timepoint)',
        'Ratio H/L (H:60 min timepoint;L:0 min. timepoint)',
        'Ratio H/L Normalized by Corresponding Protein Level (H:60 min timepoint;L:0 min. timepoint)'
    ]
]

# RM2024_names = [
#     f'RM2024_{name}' for name in [
#         'Intensity D0-1', 'Intensity D0-2', 'Intensity D0-3',
#         'Intensity D2-1', 'Intensity D2-2', 'Intensity D2-3',
#         'Intensity D4-1', 'Intensity D4-2', 'Intensity D4-3',
#         'Intensity D6-1', 'Intensity D6-2', 'Intensity D6-3',
#         'Intensity D0-1___1', 'Intensity D0-1___2', 'Intensity D0-1___3',
#         'Intensity D0-2___1', 'Intensity D0-2___2', 'Intensity D0-2___3',
#         'Intensity D0-3___1', 'Intensity D0-3___2', 'Intensity D0-3___3',
#         'Intensity D2-1___1', 'Intensity D2-1___2', 'Intensity D2-1___3',
#         'Intensity D2-2___1', 'Intensity D2-2___2', 'Intensity D2-2___3',
#         'Intensity D2-3___1', 'Intensity D2-3___2', 'Intensity D2-3___3',
#         'Intensity D4-1___1', 'Intensity D4-1___2', 'Intensity D4-1___3',
#         'Intensity D4-2___1', 'Intensity D4-2___2', 'Intensity D4-2___3',
#         'Intensity D4-3___1', 'Intensity D4-3___2', 'Intensity D4-3___3',
#         'Intensity D6-1___1', 'Intensity D6-1___2', 'Intensity D6-1___3',
#         'Intensity D6-2___1', 'Intensity D6-2___2', 'Intensity D6-2___3',
#         'Intensity D6-3___1', 'Intensity D6-3___2', 'Intensity D6-3___3'
#     ]
# ]


RN2012_names = [
    f'RN2012_H{i:02d}' for i in range(1, 7)
] + [
    f'RN2012_L{i:02d}' for i in range(10, 16)
]

RNJ2017_names = [
    f'RNJ2017_Log2_ratio_Tsup:Tstim_D{i}' for i in range(1, 4)
] + [
    f'RNJ2017_Log2_ratio_Tstim:Trest_D{i}' for i in range(2, 4)
] + [
    f'RNJ2017_Log2_ratio_Tsup:Trest_D{i}' for i in range(1, 4)
]

SAS2015_names = [
    f'SAS2015_{condition}_Exp{i}' for condition in ['AZD/DMSO', 'PLX/DMSO', 'PLX/AZD'] for i in range(1, 4)
]

SFR2015_names = [
    f'SFR2015_HL_log2_ratio_avg_{condition}' for condition in ['a', 'b', 'c']
] + [
    f'SFR2015_HL_log2_ratio_std_{condition}' for condition in ['a', 'b', 'c']
]

# SR2019_names = [
#     f'SR2019_{cell_line}_Biol_repl_{i}' for cell_line in ['BT549', 'Hs578T', 'LM2', 'MDAMB231', 'HCC1419', 'HCC1954', 
#                                                        'JIMT1', 'SKBR3', 'MCF7', 'MCF10A'] 
#     for i in range(1, 5)
# ]

TM2022_names = [
    f'TM2022_S2.{timepoint}.{sample}' 
    for timepoint in ['4h', '8h', '12h'] 
    for sample in ['A', 'B', 'C']
] + [
    f'TM2022_mock.{timepoint}.{sample}' 
    for timepoint in ['4h', '12h'] 
    for sample in ['A', 'B', 'C']
]

UK2020_names = [
    f'UK2020_HCM_{i}' for i in range(1, 7)
] + [
    f'UK2020_UA_{i}' for i in range(1, 4)
]

VOI2024_names = (
    [f'VOI2024_Veh{i} {j}' for i in range(1, 4) for j in ['mono', 'di', 'tri']] +
    [f'VOI2024_0.1uM 2 h_{i} {j}' for i in range(1, 4) for j in ['mono', 'di', 'tri']] +
    [f'VOI2024_0.1uM 24 h_{i} {j}' for i in range(1, 4) for j in ['mono', 'di', 'tri']]
)

WZ2023_names = [
    f'WZ2023_{condition}_{timepoint}_{time}min'
    for condition in ['Pe', 'Pu', 'Su']
    for timepoint in ['005', '010', '025']
    for time in ['10', '20', '30']
]

YA2020_names = [
    f'YA2020_Reporter_intensity_{condition}_No.{i}'
    for condition in ['Cancer', 'Normal']
    for i in range(1, 6)
]

YB2013_names = [
    f'YB2013_Ratio_H_L_Normalized_in_experiment_{i}' for i in range(1, 3)
]

YC2018_names = [
    f'YC2018_Log2(SKBR3-LR/SKBR3)_Rep.{rep}-{num}' 
    for rep in range(1, 3)  
    for num in range(1, 3)  
]

YP2020_names = [
    f'YP2020_Intensity {group}_{rep}_mult_log2' 
    for group in ['C10', 'C30', 'P10', 'P30']  
    for rep in ['R1', 'R2']  
]

ZCP2016_names = [f'ZCP2016_Ratio H/L Rep {i}' for i in range(1, 4)]

ZQ2022_names = [
    'ZQ2022_C1',
    'ZQ2022_C2',
    'ZQ2022_C3',
    'ZQ2022_N1',
    'ZQ2022_N2',
    'ZQ2022_N3'
]


# ZM2022_names = (
#     ['ZM2022_Intensity', 'ZM2022_Intensity___1', 'ZM2022_Intensity___2', 'ZM2022_Intensity___3'] +
#     [
#         f'ZM2022_Intensity.{cell}_{method}_{rep}'
#         for cell in ['MCF10A', 'MCF7', 'MDAMB231']
#         for method in ['IMAC_CID', 'IMAC_HCD', 'TiO2_CID', 'TiO2_HCD']
#         for rep in range(1, 4)
#     ] +
#     [
#         f'ZM2022_Intensity.{cell}_{method}_{rep}___{subrep}'
#         for cell in ['MCF10A', 'MCF7', 'MDAMB231']
#         for method in ['IMAC_CID', 'IMAC_HCD', 'TiO2_CID', 'TiO2_HCD']
#         for rep in range(1, 4)
#         for subrep in range(1, 4)
#     ]
# )

ZX2022_names = [
    f'ZX2022_Control Replicate {i}' for i in range(1, 4)
] + [
    f'ZX2022_CB treated Replicate {i}' for i in range(1, 4)
]



all_names = (# AB2020_names + 
            AH2018_names + 
            AST2020_names + AMK2021_names + AR2021_names + AS2011_names + 
            BS2014_names + CF2017_names + CK2022_names + CKC2017_names + CP2009_names + # CR2017_names + 
            CW2012_names + DB2022_names + DL2023_names + EJN2021_names + # EP2017_names + 
            FR2021_names +
            FSO2012_names + FSO2013_names + FS2019_names + GF2025_names + GRW2016_names + HH2022_names + 
            HS2024_names + HW2015_names + IMM2018_names + IMM2018_names + JAW2011_names + # JB2023_names +
            JC2019_names + JJ2018_names + JVO2010_names + JW2015_names + JW2021_names + # KBE2022_names + 
            # KDH2022_names + 
            KS2014_names + 
            # KTGR2014_names + 
            LATM2013_names + LG2023_names + LY2019_names +
            MM2018_names + MV2014_names + NJH2015_names + NJH2024_names + NJH2025_names + PG2020_names + 
            RAM2015_names + RB2022_names + RBK2021_names + RCJD2014_names + RKK2015_names + RM2009_names + 
            # RM2024_names + 
            RN2012_names + RNJ2017_names + SAS2015_names + SFR2015_names + # SR2019_names + 
            TM2022_names + UK2020_names + VOI2024_names + WZ2023_names + 
            YA2020_names + YB2013_names + 
            YC2018_names + YP2020_names + ZCP2016_names + # ZM2022_names + 
            ZQ2022_names + ZX2022_names)


# ----------------- #
# CREATE DICTIONARY OF DATAFRAMES
# ----------------- #

dataset_dict = {}

# add each dataframe from list to dictionary (with index)
for i, dataset_name in enumerate(all_names):
    dataset_dict[dataset_name] = dataset_list[i]

print(dataset_dict[f'AH2018_TMTMS2_1_C'].head())  


# ----------------- #
# CREATE MATRIX HEADER
# ----------------- #

norm_matrix = normalising.MinMax_normalize_and_merge(dataset_dict, MinMaxScaler())

norm_matrix.columns.name = ''

# get unique phosphosites
uniq_phos = norm_matrix.columns.drop('DatasetName')

if len(uniq_phos) == len(norm_matrix.columns) - 1:
    print('Passed! Number of unique phosphosites matches expected number of columns')
else:
    print('Failed! Number of unique phosphosites does not match expected number of columns')

matrix_cols = pd.DataFrame(columns = uniq_phos)

matrix_cols.to_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/normalised-matrix-header.csv', index = False)

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

minmax_matrix.to_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NormalisedMatrix.csv-NBA', index = False)
print(f'Normalised matrix saved successfully!', minmax_matrix)
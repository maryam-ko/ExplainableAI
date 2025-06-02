import pandas as pd
import numpy as np
from scipy.stats import zscore

csv = pd.read_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NormalisedMatrix-NBA.csv')

csv['MedianCoeff'] = csv.median(axis=1, numeric_only=True)

csv['MedianCoeff_z'] = zscore(csv['MedianCoeff'])

# Filter out datasets where |z-score| >= 3 
filtered_csv = csv[np.abs(csv['MedianCoeff_z']) < 3]

filtered_csv = filtered_csv.drop(columns=['MedianCoeff', 'MedianCoeff_z'])

filtered_csv.to_csv(
    '/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NormalisedMatrix-Zscore.csv',
    index=False
)

print(f"Original datasets: {csv.shape[0]}")
print(f"Datasets after filtering: {filtered_csv.shape[0]}")
print("Z-score filtering complete.")

#!/bin/python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load normalized phosphosite data (samples x phosphosites)
df = pd.read_csv('/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NormalisedMatrix-Zscore.csv', index_col=0)
print(df)
# Clean NaN and infinite values
df = df.replace([np.nan, np.inf, -np.inf], 0)
print(df)
# Extract protein names from phosphosite column names (assumes 'PROTEIN_Site' format)
protein_names = df.columns.str.split('_').str[0]
print(protein_names)
# Create a DataFrame mapping phosphosites to proteins
phosphosite_to_protein = pd.DataFrame({'phosphosite': df.columns, 'protein': protein_names})
print(phosphosite_to_protein)
# Group columns by protein: For each protein, average all phosphosite columns
protein_level_df = pd.DataFrame(index=df.index)  
print(protein_level_df)

for protein in phosphosite_to_protein['protein'].unique():
    cols = phosphosite_to_protein.loc[phosphosite_to_protein['protein'] == protein, 'phosphosite']
    protein_level_df[protein] = df[cols].mean(axis=1)

print(f"Protein-level DataFrame shape: {protein_level_df.shape}")
print(protein_level_df.head())
print(protein_level_df)

# Compute Spearman correlation matrix between proteins
corr_matrix = protein_level_df.corr(method='spearman').fillna(0)
print(corr_matrix)

# Plot the heatmap
plt.figure(figsize=(12, 10))
ax = sns.heatmap(corr_matrix, cmap='coolwarm', center=0, xticklabels=False, yticklabels=False)
ax.set_xlabel('')
ax.set_ylabel('')
plt.tight_layout()
plt.savefig('protein_correlation_heatmap1.png', dpi=300)
plt.show()

# Plot the clustered heatmap with dendrograms
sns.clustermap(
    corr_matrix,
    cmap='coolwarm',
    center=0,
    method='average',     # Clustering method
    metric='euclidean',   
    figsize=(12, 12),
    xticklabels=False,   
    yticklabels=False
).savefig('protein_correlation_clustermap_with_dendrograms1.png', dpi=300)
plt.show()

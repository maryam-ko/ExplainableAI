#!/bin/python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load your matrix
df = pd.read_csv("/data/home/bt24990/ExplainableAI/04_clustering/normalised_matrix.csv", index_col=0)

# Columns are phosphosites like 'AKT1_T308'
phosphosites = df.columns.tolist()
num_phosphosites = len(phosphosites)
print(f"Number of phosphosites: {num_phosphosites}")

# Extract protein names from phosphosite labels
proteins = [p.split('_')[0] for p in phosphosites]
unique_proteins = set(proteins)
num_proteins = len(unique_proteins)
print(f"Number of unique proteins: {num_proteins}")

# Plot bar chart with new axis labels
plt.figure(figsize=(6, 5))
plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)
bars = plt.bar(['Phosphosites', 'Proteins'], [num_phosphosites, num_proteins], color=['slateblue', 'slateblue'], width=0.3, zorder=3)

# Add value labels on top of bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2,
             height + 5,
             f'{int(height)}',
             ha='center', va='bottom',
             fontsize=12, fontweight='bold')


# Axis labels
plt.xlabel('Threshold', fontsize=14, fontweight='bold', labelpad=10)
plt.ylabel('Number of features', fontsize=14, fontweight='bold', labelpad=10)

plt.tight_layout(pad=2.0)
plt.show()
plt.savefig("bar_plot.png", dpi=300, bbox_inches='tight', pad_inches=0)


#!/bin/python

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

paths = {
    '0': "/data/home/bt24990/ExplainableAI/04_clustering/interim_data/clustered_matrix.csv",
    '50': "/data/home/bt24990/ExplainableAI/04_clustering/clustered_matrix_min50vals.csv",
    '100': "/data/home/bt24990/ExplainableAI/04_clustering/clustered_matrix_min100vals.csv",
    '150': "/data/home/bt24990/ExplainableAI/04_clustering/clustered_matrix_min150vals.csv",
    '200': "/data/home/bt24990/ExplainableAI/04_clustering/clustered_matrix_min200vals.csv"
}

counts = {}
for label, path in paths.items():
    df = pd.read_csv(path, index_col=0)
    counts[label] = len(df.columns)
    print(f"Number of phosphosites in {label} matrix: {counts[label]}")

# Plot bar chart with new axis labels
plt.figure(figsize=(7, 5))
plt.grid(axis='y', linestyle='--', alpha=0.7, zorder=0)

bars = plt.bar(counts.keys(), counts.values(), width=0.4, zorder=3, color='blue')

# Add value labels
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width() / 2,
             height + 5,
             f'{int(height)}',
             ha='center', va='bottom',
             fontsize=12, fontweight='bold')

# Axis labels
plt.xlabel('Thresholds', fontsize=14, fontweight='bold', labelpad=10)
plt.ylabel('Number of features', fontsize=14, fontweight='bold', labelpad=10)

plt.tight_layout(pad=2.0)
plt.savefig("/data/home/bt24990/ExplainableAI/08_results/linear_regression/plots/barplot_features.png",
            dpi=300, bbox_inches='tight', pad_inches=0)
plt.show()
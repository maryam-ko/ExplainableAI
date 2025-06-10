#!/bin/python

import os
import pandas as pd

input_dir = "/data/home/bt24990/ExplainableAI/04_clustering/interim_data"
output_file = "/data/home/bt24990/ExplainableAI/04_clustering/combined_clustered_matrix.csv"

dfs = []
for file in os.listdir(input_dir):
    if file.endswith("_clustered_matrix.csv"):
        df = pd.read_csv(os.path.join(input_dir, file), index_col=0)
        dfs.append(df)

# Concatenate side by side (same rows, different columns)
combined = pd.concat(dfs, axis=1)

# Sort columns alphabetically
combined = combined.sort_index(axis=1)

# Save to CSV
combined.to_csv(output_file, index=True)
print(f"Combined matrix saved to {output_file}")
print(f"Combined matrix shape: {combined.shape}  # (rows, columns)")

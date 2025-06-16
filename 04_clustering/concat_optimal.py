#!/bin/python

import pandas as pd
import glob
import os

# Directory containing the files
interim_data_dir = '/data/home/bt24990/ExplainableAI/04_clustering/interim_data'

# Find all *_optimal_clusters.csv files
csv_files = glob.glob(os.path.join(interim_data_dir, '*_optimal_clusters.csv'))

# Read and concatenate
dfs = []
for file in csv_files:
    df = pd.read_csv(file)
    dfs.append(df)

# Concatenate along rows
all_clusters = pd.concat(dfs, axis=0, ignore_index=True)

# # Fill empty cells with 0.0 (important: do this before reordering columns)
# all_clusters.fillna(0.0, inplace=True)

# Reorder columns: Prefix, ModeClusters, then all random seed columns
cols = ['Prefix', 'ModeClusters'] + [col for col in all_clusters.columns if col not in ['Prefix', 'ModeClusters']]
all_clusters = all_clusters[cols]

# Save to CSV
all_clusters.to_csv('/data/home/bt24990/ExplainableAI/04_clustering/all_optimal_clusters.csv', index=False)

print(f"Concatenated {len(csv_files)} files into all_optimal_clusters.csv")
print(f"Total rows in concatenated DataFrame: {len(all_clusters)}")
print(f"Columns in concatenated DataFrame: {all_clusters.columns.tolist()}")
print(all_clusters.head())

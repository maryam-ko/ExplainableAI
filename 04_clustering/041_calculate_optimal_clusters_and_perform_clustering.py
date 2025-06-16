#!/bin/python

import time
start_time = time.time()
import pandas as pd
import numpy as np
import os
import sys
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from multiprocessing import Pool

grandparent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(grandparent_dir)

from funcs import generalfuncs

import argparse

# ----------------- #

# def process_file(file):
#     split_dir = '/data/home/bt24990/ExplainableAI/04_clustering/split_matrices'
#     output_dir = '/data/home/bt24990/ExplainableAI/04_clustering/interim_data'
#     protein_name = os.path.basename(file).replace('_matrix.csv', '')
#     output_file = os.path.join(output_dir, f"{protein_name}_optimal_clusters.csv")

#     if os.path.exists(output_file):
#         print(f"Skipping {protein_name} — output already exists.")
#         return

#     input_file = os.path.join(split_dir, file)
#     print(f"Processing: {protein_name}")

#     matrix = pd.read_csv(input_file, index_col=0)
#     matrix = generalfuncs.set_dataset_name_as_index(matrix)

#     print("Calculating the optimal number of clusters over 100 random seeds...")
#     try:
#         optimal_clusters = generalfuncs.calculate_optimal_clusters(
#             matrix,
#             filename=f"{protein_name}_optimal_clusters"
#         )
#         print(f"Finished processing {protein_name}\n")
#     except ValueError as e:
#         print(f"Skipping {protein_name}: {e}")

def process_file(file):
    split_dir = '/data/home/bt24990/ExplainableAI/04_clustering/split_matrices'
    interim_dir = '/data/home/bt24990/ExplainableAI/04_clustering/interim_data'
    protein_name = os.path.basename(file).replace('_matrix.csv', '')

    optimal_cluster_file = os.path.join(interim_dir, f"{protein_name}_optimal_clusters.csv")

    # Skip if optimal cluster file doesn't exist
    if not os.path.exists(optimal_cluster_file):
        print(f"Skipping {protein_name} — no optimal_clusters file.")
        return

    print(f"Processing: {protein_name}")
    input_file = os.path.join(split_dir, file)
    matrix = pd.read_csv(input_file, index_col=0)
    matrix = generalfuncs.set_dataset_name_as_index(matrix)

    optimal_clusters = pd.read_csv(optimal_cluster_file, index_col=0)

    try:
        generalfuncs.create_clustered_matrix_from_normalised_matrix(
            matrix,
            optimal_clusters,
            filename=f"{protein_name}_clustered_matrix"
        )
        print(f"Finished clustering {protein_name}\n")
    except Exception as e:
        print(f"Failed to cluster {protein_name}: {e}")

if __name__ == "__main__":
    split_dir = '/data/home/bt24990/ExplainableAI/04_clustering/split_matrices'
    files = [f for f in os.listdir(split_dir) if f.endswith('_matrix.csv')]

    n_processes = 8
    with Pool(processes=n_processes) as pool:
        pool.map(process_file, files)

    print(f'Execution time: {time.time() - start_time:.2f} seconds, {(time.time() - start_time)/60:.2f} minutes, {(time.time() - start_time)/3600:.2f} hours.')
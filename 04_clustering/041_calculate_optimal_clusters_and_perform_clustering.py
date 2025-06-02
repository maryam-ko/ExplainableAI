#!/bin/python

import time
start_time = time.time()
import pandas as pd
import numpy as np
import os
import sys
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

grandparent_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(grandparent_dir)

from funcs import generalfuncs

# add argparse function to allow for command line arguments and specifically, the base directory

# ----------------- #

if __name__ == "__main__":
    
    print("Loading normalised matrix...")
    # import the phosphosite-resolved and normalised matrix
    matrix = pd.read_csv('/data/home/bt24990/ExplainableAI/04_clustering/normalised_matrix.csv', header=0)
    matrix = generalfuncs.set_dataset_name_as_index(matrix)
    print(matrix.head())

    print("Calculating the optimal number of clusters per protein over 100 random seeds...")    
    optimal_clusters = generalfuncs.calculate_optimal_clusters(matrix, 
                                                               filename="optimal_clusters")
    
    print("Performing clustering to group phosphosites...")
    clustered_matrix = generalfuncs.create_clustered_matrix_from_normalised_matrix(matrix, 
                                                                                   optimal_clusters,
                                                                                   filename="clustered_matrix")
    print("Clustering complete. Clustered matrix has been saved.")


    print(f'Execution time: {time.time() - start_time:.2f} seconds, {(time.time() - start_time)/60:.2f} minutes, {(time.time() - start_time)/3600:.2f} hours.')
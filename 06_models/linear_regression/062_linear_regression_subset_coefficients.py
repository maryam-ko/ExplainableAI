#!/bin/python

import time
start_time = time.time()
import pandas as pd
import os
import sys
import argparse

# ----------------- #

def parse_arguments():
    """Parse command line arguments for the script."""

    parser = argparse.ArgumentParser(description='Evaluate protein interaction predictions against Biogrid reference.')
    parser.add_argument('--base_dir', type=str, default='/data/home/bt24990/ExplainableAI', help='Base directory for the project')
    parser.add_argument('--threshold', type=int, default=200, help='Threshold to compute')
    parser.add_argument('--network_name', type=str, default='MAPKERK', help='Network name to compute')

    return parser.parse_args()

def subset_lr_coefficients(base_dir, threshold, network_name):
    """Subset LR coefficients to only include specified proteins."""

    print(f"Subsetting LR coefficients for {network_name} at threshold = {threshold}...")
    df = pd.read_csv(f"{base_dir}/08_results/linear_regression/coefficients/linear_regression_cv_coefficients_min{threshold}vals.csv", header=0)
    df.loc[:, 'Feature'] = df['Feature'].apply(lambda x: x.split("_")[0])
    df.loc[:, 'TargetFeature'] = df['TargetFeature'].apply(lambda x: x.split("_")[0])
    df_subset = df[df.apply(lambda row: row['Feature'] in prots or row['TargetFeature'] in prots, axis=1)]
    df_subset.to_csv(f"{base_dir}/08_results/linear_regression/coefficients/linear_regression_nested_cv_{network_name}_coefficients_protein_level_min{threshold}vals.csv", index=False)

# ----------------- #

if __name__ == '__main__':
    
    args = parse_arguments()

    if args.network_name == "MAPKERK":
        prots = (
            'ARAF', 'ATF2', 'BRAF', 'EGFR', 'ELK1', 
            'ERBB2','ETS1', 'GAB1', 'FOS', 'FRS2', 
            'IRS1','MAPK1','MAPK3', 'MAP2K1','MAP2K2', 
            'MAPKAPK2','MET', 'MEF2A', 'MEF2C', 'MITF',
            'PDGFRA', 'RAF1', 'RPS6KA1','RPS6KA3', 'SOS1', 
            'SOS2', 'SHC1', 'SPRED1', 'SPRED2'
        )

    subset_lr_coefficients(
        base_dir=args.base_dir, 
        threshold=args.threshold,
        network_name=args.network_name
    )

    print(f'Execution time: {time.time() - start_time:.2f} seconds, {(time.time() - start_time)/60:.2f} minutes, {(time.time() - start_time)/3600:.2f} hours.')
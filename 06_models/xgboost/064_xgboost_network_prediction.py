#!/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyvis.network import Network
import time
start_time = time.time()
import os
import sys
import argparse

grandparent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.append(grandparent_dir)
from funcs import plots

network_name = "MAPKERK"
prots = ('ARAF', 'ATF2', 'BRAF', 'EGFR', 'ELK1', 
            'ERBB2','ETS1', 'GAB1', 'FOS', 'FRS2', 
            'IRS1','MAPK1','MAPK3', 'MAP2K1','MAP2K2', 
            'MAPKAPK2','MET', 'MEF2A', 'MEF2C', 'MITF',
            'PDGFRA', 'RAF1', 'RPS6KA1','RPS6KA3', 'SOS1', 
            'SOS2', 'SHC1', 'SPRED1', 'SPRED2')

# ----------------- #

def parse_arguments():
    """Parse command line arguments for the script."""

    parser = argparse.ArgumentParser(description='Evaluate protein interaction predictions against Biogrid reference.')
    parser.add_argument('--base_dir', type=str, default='/data/home/bt24990/ExplainableAI', help='Base directory for the project')
    parser.add_argument('--threshold', type=int, default=150, help='Threshold to compute')

    return parser.parse_args()


# ----------------- #

if __name__ == "__main__":

    args = parse_arguments()

    # print('Loading linear regression predicted interactions...')
    # lr_df = pd.read_csv(f"{args.base_dir}/08_results/linear_regression/coefficients/linear_regression_nested_cv_coefficients_min{args.threshold}vals.csv")

    # # protein level
    # lr_df['TargetFeature'] = lr_df['TargetFeature'].str.split('_').str[0]
    # lr_df['PredictiveFeature'] = lr_df['PredictiveFeature'].str.split('_').str[0]
    # lr_df = lr_df[['PredictiveFeature', 'TargetFeature', 'MedianCoeff']]
    # print('Protein level linear regression predictions:', lr_df)

    # plots.plot_LR_predicted_network(
    #     args.base_dir,
    #     lr_df, 
    #     strong_boundary=0.57, 
    #     medium_boundary=0.55,
    #     weak_boundary=0.53, 
    #     selected_prots=prots, 
    #     save_as_filename=f"linear_regression_nested_cv_predicted_network_min{args.threshold}vals.html")
    
    # ----------------- #
    
    print('Selecting MAPERK proteins...')
    subset_filename = f"xgboost_master_shap_file_cluster_level_min{args.threshold}vals.csv"
    df = pd.read_csv(f'{args.base_dir}/08_results/xgboost/nested_cv_master_shaps/{subset_filename}')
    df['TargetFeature'] = df['TargetFeature'].str.split('_').str[0]
    df['PredictiveFeature'] = df['Feature'].str.split('_').str[0]
    print('Protein level xgboost predictions:', df)
    subset_df = df[df.apply(lambda row: row['Feature'] in prots or row['TargetFeature'] in prots, axis=1)]
    print('Subset xgboost predictions:', subset_df)

    plots.plot_CNN_or_XGB_predicted_network(
        base_dir=args.base_dir,
        model='xgboost',
        df=subset_df, 
        strong_percentile=1, 
        medium_percentile=3,
        weak_percentile=5, 
        selected_prots=prots,
        save_as_filename=f"xgboost_nested_cv_{network_name}_predicted_network_min{args.threshold}vals.html")
    
    print(f'Execution time: {time.time() - start_time:.2f} seconds, {(time.time() - start_time)/60:.2f} minutes, {(time.time() - start_time)/3600:.2f} hours.')



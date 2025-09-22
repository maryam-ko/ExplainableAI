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
    parser.add_argument('--threshold', type=int, default=200, help='Threshold to compute')

    return parser.parse_args()

def plot_XGB_predicted_network(base_dir, df, strong_percentile, medium_percentile, 
                              weak_percentile, selected_prots, save_as_filename):
    """
    Plots network of selected proteins using XGBoost SHAP values.

    Inputs:
    df <dataframe>: DataFrame with columns ['PredictiveFeature', 'TargetFeature', 'SHAPValue', 'SHAP*R2']
    strong_percentile, medium_percentile, weak_percentile <int>: Percentile cutoffs for edge strength
    selected_prots <tuple>: Proteins to highlight
    save_as_filename <str>: Output filename
    """
    abs_shap = df['SHAP*R2']
    strong_boundary = np.percentile(abs_shap, 100 - strong_percentile)
    medium_boundary = np.percentile(abs_shap, 100 - medium_percentile)
    weak_boundary = np.percentile(abs_shap, 100 - weak_percentile)

    nt = Network('1000px', '1000px', directed=True)
    for row in df.itertuples(index=False):
        PredictiveFeature = row.PredictiveFeature
        TargetFeature = row.TargetFeature
        SHAPValue = row._2  # or row.SHAP_xR2 if you rename the column

        # Set default colors
        feature_color = '#a7b5e0'
        targetfeat_color = '#a7b5e0'

        # Highlight selected proteins
        if PredictiveFeature in selected_prots:
            feature_color = '#d16002'
        if TargetFeature in selected_prots:
            targetfeat_color = '#d16002'

        strong_positive = '#36441d'
        medium_positive = '#708d3e'
        weak_positive = '#d4e2bd'
        strong_negative = '#873c07'
        medium_negative = '#ca5a0b'
        weak_negative = '#f8bb8f'
        node_fontsize = 25

        # Feedback loops
        if PredictiveFeature == TargetFeature:
            if SHAPValue <= -strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=strong_negative, value=2)
            elif SHAPValue <= -medium_boundary and SHAPValue > -strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=medium_negative, value=2)
            elif SHAPValue <= -weak_boundary and SHAPValue > -medium_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=weak_negative, value=2)
            elif SHAPValue >= strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=strong_positive, value=2)
            elif SHAPValue >= medium_boundary and SHAPValue < strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=medium_positive, value=2)
            elif SHAPValue >= weak_boundary and SHAPValue < medium_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=weak_positive, value=2)
        # All other edges
        else:
            if SHAPValue <= -strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=strong_negative, value=2)
            elif SHAPValue <= -medium_boundary and SHAPValue > -strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=medium_negative, value=2)
            elif SHAPValue <= -weak_boundary and SHAPValue > -medium_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=weak_negative, value=2)
            elif SHAPValue >= strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=strong_positive, value=2)
            elif SHAPValue >= medium_boundary and SHAPValue < strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=medium_positive, value=2)
            elif SHAPValue >= weak_boundary and SHAPValue < medium_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=weak_positive, value=2)

    nt.save_graph(f'{base_dir}/08_results/xgboost/predicted_networks/{save_as_filename}')
    print('Predicted XGBoost network graph saved!')

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
    subset_filename = f"xgboost_master_shap_files_cluster_level_min{args.threshold}vals_shap_r2_short.csv"
    df = pd.read_csv(f'{args.base_dir}/06_models/xgboost/master_shaps_files/{subset_filename}')
    df['TargetFeature'] = df['TargetFeature'].str.split('_').str[0]
    df['PredictiveFeature'] = df['PredictiveFeature'].str.split('_').str[0]
    print('Protein level xgboost predictions:', df)
    subset_df = df[df.apply(lambda row: row['PredictiveFeature'] in prots or row['TargetFeature'] in prots, axis=1)]
    subset_df = subset_df[['PredictiveFeature', 'TargetFeature', 'SHAP*R2']]
    print('Subset xgboost predictions:', subset_df)

    plot_XGB_predicted_network(
        args.base_dir,
        subset_df, 
        strong_percentile=1, 
        medium_percentile=3,
        weak_percentile=5, 
        selected_prots=prots, 
        save_as_filename=f"xgboost_nested_cv_{network_name}_predicted_network_min{args.threshold}vals.html")

    print(f'Execution time: {time.time() - start_time:.2f} seconds, {(time.time() - start_time)/60:.2f} minutes, {(time.time() - start_time)/3600:.2f} hours.')



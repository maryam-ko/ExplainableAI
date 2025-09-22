#!/bin/python

# Run from command line
# python3 omnipath_evaluation.py --base_dir /path/to/project --network_name my_network --model_types linear_regression xgboost --thresholds 50 100
# or just
# python3 omnipath_evaluation.py 
# for default parameters

import time
start_time = time.time()
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve, roc_curve, auc, roc_auc_score
import os
import argparse


# ----------------- #

def parse_arguments():
    """Parse command line arguments for the script."""
    parser = argparse.ArgumentParser(description='Evaluate protein interaction predictions against OmniPath reference.')
    parser.add_argument('--base_dir', type=str, default='/data/home/bt24990/ExplainableAI', help='Base directory for the project')
    parser.add_argument('--network_name', type=str, default=None, help='Network name if analyzing specific network')
    parser.add_argument('--model_types', nargs='+', default=['linear_regression'], 
                        help='Model types to evaluate')
    parser.add_argument('--thresholds', nargs='+', type=int, default=[50, 100, 150, 200], 
                        help='Thresholds to evaluate')
    return parser.parse_args()

def get_file_path(base_dir, network_name, model_type, threshold, file_type, int_type):
    # """Construct file paths based on parameters."""
    # if file_type == 'input_preds':
    #     if model_type == 'linear_regression':
    #         return f'{base_dir}/08_results/linear_regression/coefficients/linear_regression_min{threshold}vals_coeff_r2.csv'
    #     elif model_type == 'xgboost':
    #         return f'{base_dir}/06_models/{model_type}/master_shaps_files/{model_type}_master_shap_files_cluster_level_min{threshold}vals_shap_r2_short.csv'
    # elif file_type == 'csv_output':
    #     filename = f'omnipath_aucs_{int_type}'
    #     if network_name is not None:
    #         filename = f'{filename}_{network_name}'
    #     return f'{base_dir}/08_results/omnipath_results/{filename}.csv'
    # else:
    #     # For plot outputs
    #     filename = f'omnipath_aucs_{int_type}'
    #     if network_name is not None:
    #         filename = f'{filename}_{network_name}'
    #     return f'{base_dir}/08_results/omnipath_results/{filename}.png'
    """Construct file paths based on parameters."""
    if file_type == 'input_preds':
        if model_type == 'linear_regression':
            if network_name is None:
                return f'{base_dir}/08_results/linear_regression/coefficients/linear_regression_cv_coefficients_min{threshold}vals.csv'
            else:
                return f'{base_dir}/08_results/linear_regression/coefficients/linear_regression_nested_cv_{network_name}_coefficients_protein_level_min{threshold}vals.csv'
        else:
            if network_name is None:
                return f'{base_dir}/06_models/{model_type}/master_shaps_files/{model_type}_master_shap_file_cluster_level_min{threshold}vals_shap_r2.csv'
            else:
                return f'{base_dir}/06_models/{model_type}/master_shaps_files/{model_type}_master_shap_file_protein_level_min{threshold}vals.csv'
    elif file_type == 'csv_output':
        filename = f'omnipath_aucs_{int_type}'
        if network_name is not None:
            filename = f'{filename}_{network_name}'
        return f'{base_dir}/08_results/omnipath_results/{filename}.csv'
    else:
        # For plot outputs
        filename = f'omnipath_aucs_{int_type}'
        if network_name is not None:
            filename = f'{filename}_{network_name}'
        return f'{base_dir}/08_results/omnipath_results/{filename}.png'
    
    
def load_omnipath_reference(base_dir):
    df = pd.read_csv(f"{base_dir}/07_evaluation/omnipath_evaluation/omnipath_human_interactions.csv", sep='\t')
    df = df[['source_gene_name', 'target_gene_name']]
    df.columns = ['PredictiveFeature', 'TargetFeature']
    return df

def load_biogrid_reference(base_dir):
    df = pd.read_csv(f"{base_dir}/01_input_data/BIOGRID-ORGANISM-Homo_sapiens-4.4.246.tab3.txt", sep="\t", header=0, low_memory=False)
    df = df[['Official Symbol Interactor A', 'Official Symbol Interactor B']]
    df2 = df.rename(columns={'Official Symbol Interactor A': 'Official Symbol Interactor B', 'Official Symbol Interactor B': 'Official Symbol Interactor A'})
    df = pd.concat([df, df2], axis=0)
    df.columns = ['PredictiveFeature', 'TargetFeature']
    df.drop_duplicates(inplace=True)
    return df

def load_prediction_data(base_dir, network_name, model_type, threshold):
    """Load and format the predicted interactions data."""
    file_path = get_file_path(base_dir=base_dir, 
                              network_name=network_name, 
                              model_type=model_type, 
                              threshold=threshold,
                              file_type='input_preds', 
                              int_type=None)
    print(f"Loading prediction data from: {file_path}")
    df = pd.read_csv(file_path, header=0)
    
    if model_type == 'linear_regression':
        df['TargetFeature'] = df['TargetFeature'].str.split('_').str[0]
        df['PredictiveFeature'] = df['PredictiveFeature'].str.split('_').str[0]
        df = df[['PredictiveFeature', 'TargetFeature', 'MedianCoeff']]
        value_col = 'MedianCoeff'

    # if model_type == 'linear_regression':
    #     df['TargetFeature'] = df['TargetFeature'].str.split('_').str[0]
    #     df['PredictiveFeature'] = df['PredictiveFeature'].str.split('_').str[0]
    #     df = df[['PredictiveFeature', 'TargetFeature', 'Coeff*R2']]
    #     value_col = 'Coeff*R2'
    else:
        if 'SHAPValue' in df.columns:
            df['PredictiveFeature'] = df['PredictiveFeature'].apply(lambda x: ' '.join(sorted(set([token.split('_')[0] for token in x.split(' ')]))))
            df['TargetFeature'] = df['TargetFeature'].apply(lambda x: ' '.join(sorted(set([token.split('_')[0] for token in x.split(' ')]))))
            df = df[['PredictiveFeature', 'TargetFeature', 'SHAPValue']]
            value_col = 'SHAPValue'
    # else:
    #     if 'SHAP*R2' in df.columns:
    #         df['PredictiveFeature'] = df['PredictiveFeature'].apply(lambda x: ' '.join(sorted(set([token.split('_')[0] for token in x.split(' ')]))))
    #         df['TargetFeature'] = df['TargetFeature'].apply(lambda x: ' '.join(sorted(set([token.split('_')[0] for token in x.split(' ')]))))
    #         df = df[['PredictiveFeature', 'TargetFeature', 'SHAP*R2']]
    #         value_col = 'SHAP*R2'
        # else:
        #     raise KeyError("No SHAPValue or SHAP*R2 column found in XGBoost predictions.")
    
    # Remove duplicates
    df.dropna(subset=[value_col], inplace=True)
    df.drop_duplicates(subset=['PredictiveFeature', 'TargetFeature'], inplace=True)
    print('Predicted dataframe:', df)
    return df, value_col

def filter_for_common_proteins(predictions, reference):
    """Retain only proteins present in both datasets."""
    # Standardize protein names
    for df in [predictions, reference]:
        df['PredictiveFeature'] = df['PredictiveFeature'].str.upper().str.strip()
        df['TargetFeature'] = df['TargetFeature'].str.upper().str.strip()
    
    # Get common proteins
    pred_prots = set(predictions['PredictiveFeature']) | set(predictions['TargetFeature'])
    ref_prots = set(reference['PredictiveFeature']) | set(reference['TargetFeature'])
    common_prots = pred_prots & ref_prots
    print(f'Number of common proteins: {len(common_prots)}')
    
    # Filter both datasets
    mask_preds = (predictions['PredictiveFeature'].isin(common_prots)) & (predictions['TargetFeature'].isin(common_prots))
    mask_ref = (reference['PredictiveFeature'].isin(common_prots)) & (reference['TargetFeature'].isin(common_prots))
    
    return predictions[mask_preds].copy(), reference[mask_ref].copy()

def format_df_per_model_and_threshold(base_dir, network_name, model_types, thresholds):
    """Main function to run the evaluation process."""
    
    # Ensure output directory exists
    os.makedirs(f'{base_dir}/08_results/omnipath_results', exist_ok=True)
    
    print("Loading and processing prediction data...")
    predictions_dict = {model: {} for model in model_types}
    value_cols = {model: None for model in model_types}
    reference_filtered_dict = {model: {} for model in model_types}

    print('Creating dictionaries of predictions and reference dataframes...')
    for model in model_types:
        # Load the correct reference for each model
        if model == 'linear_regression':
            reference_data = load_biogrid_reference(base_dir)
        else:
            reference_data = load_omnipath_reference(base_dir)
        for threshold in thresholds:
            # Load predictions
            predictions, value_col = load_prediction_data(base_dir, network_name, model, threshold)
            value_cols[model] = value_col
            
            # Filter for common proteins
            filtered_preds, filtered_ref = filter_for_common_proteins(predictions, reference_data)
            predictions_dict[model][threshold] = filtered_preds
            reference_filtered_dict[model][threshold] = filtered_ref
    
    return predictions_dict, reference_filtered_dict, value_cols

def calculate_aucs(preds_dict, ref_dict, value_cols, base_dir, network_name, model_types, thresholds, int_type):
    """Calculate PR-AUC and ROC-AUC scores."""
    auc_results = []
    for model in model_types:
        model_results = []
        for threshold in thresholds:
            preds_dict[model][threshold] = preds_dict[model][threshold].drop_duplicates(subset=['PredictiveFeature', 'TargetFeature'])
            ref_dict[model][threshold] = ref_dict[model][threshold].drop_duplicates(subset=['PredictiveFeature', 'TargetFeature'])

            y_scores = preds_dict[model][threshold][value_cols[model]].tolist()

            merged_df = pd.merge(preds_dict[model][threshold], ref_dict[model][threshold], on=['PredictiveFeature', 'TargetFeature'],
                                how='left', indicator=True)
            y_true = (merged_df['_merge'] == 'both').astype(int).tolist()

            precision, recall, _ = precision_recall_curve(y_true, y_scores)
            pr_auc = auc(recall, precision)
            roc_auc = roc_auc_score(y_true, y_scores)
            print(f"Positive vs negative examples ({model}, {threshold}): {sum(y_true)} vs {len(y_true) - sum(y_true)}")
            
            model_results.append({
                'model_type': model,
                'threshold': threshold,
                'pr_auc': pr_auc,
                'roc_auc': roc_auc,
                'positive_examples': sum(y_true),
                'negative_examples': len(y_true) - sum(y_true)
            })
        auc_results.extend(model_results)

    auc_df = pd.DataFrame(auc_results)
    output_result_path = get_file_path(
        base_dir=base_dir,
        network_name=network_name, 
        model_type=None, 
        threshold=None, 
        file_type='csv_output',
        int_type=int_type
    )
    auc_df.to_csv(output_result_path, index=False)
    print("AUC results:", auc_df)

    return auc_df


# ----------------- #

if __name__ == '__main__':

    args = parse_arguments()

    preds_dict, ref_dict, value_cols = format_df_per_model_and_threshold(
        base_dir=args.base_dir,
        network_name=args.network_name, 
        model_types=args.model_types, 
        thresholds=args.thresholds
    )

    # --- Loop over all thresholds ---
    for threshold in args.thresholds:

        # Linear Regression: OmniPath vs BioGRID
        if 'linear_regression' in preds_dict:
            preds_lr = preds_dict['linear_regression'][threshold]

            # Against BioGRID
            ref_biogrid = load_biogrid_reference(args.base_dir)
            preds_lr_bio, ref_bio = filter_for_common_proteins(preds_lr, ref_biogrid)
            merged_lr_bio = pd.merge(
                preds_lr_bio,
                ref_bio,
                on=['PredictiveFeature', 'TargetFeature'],
                how='left',
                indicator=True
            )
            y_true_lr_bio = (merged_lr_bio['_merge'] == 'both').astype(int).tolist()
            y_scores_lr_bio = merged_lr_bio[value_cols['linear_regression']].tolist()
            fpr_lr_bio, tpr_lr_bio, _ = roc_curve(y_true_lr_bio, y_scores_lr_bio)
            roc_auc_lr_bio = auc(fpr_lr_bio, tpr_lr_bio)

            # Against OmniPath
            ref_omnipath = load_omnipath_reference(args.base_dir)
            preds_lr_omni, ref_omni = filter_for_common_proteins(preds_lr, ref_omnipath)
            merged_lr_omni = pd.merge(
                preds_lr_omni,
                ref_omni,
                on=['PredictiveFeature', 'TargetFeature'],
                how='left',
                indicator=True
            )
            y_true_lr_omni = (merged_lr_omni['_merge'] == 'both').astype(int).tolist()
            y_scores_lr_omni = merged_lr_omni[value_cols['linear_regression']].tolist()
            fpr_lr_omni, tpr_lr_omni, _ = roc_curve(y_true_lr_omni, y_scores_lr_omni)
            roc_auc_lr_omni = auc(fpr_lr_omni, tpr_lr_omni)

            # Plot both ROC curves together for Linear Regression
            plt.figure(figsize=(7, 6))
            plt.plot(fpr_lr_omni, tpr_lr_omni, color='orange', lw=2, label=f'LR OmniPath ROC (AUC = {roc_auc_lr_omni:.2f})')
            plt.plot(fpr_lr_bio, tpr_lr_bio, color='blue', lw=2, label=f'LR BioGRID ROC (AUC = {roc_auc_lr_bio:.2f})')
            plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate', fontsize=14)
            plt.ylabel('True Positive Rate', fontsize=14)
            plt.title(f'Linear Regression ROC: OmniPath vs BioGRID\n(Threshold {threshold}, CoeffValues)', fontsize=15)
            # plt.title(f'Linear Regression ROC: OmniPath vs BioGRID\n(Threshold {threshold}, Coeff*R2)', fontsize=15)
            plt.legend(loc="lower right")
            plt.grid(True, linestyle='--', alpha=0.3)
            plt.tight_layout()
            plt.savefig(f"{args.base_dir}/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold{threshold}.png")
            # plt.savefig(f"{args.base_dir}/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold{threshold}_coeff_r2.png")
            plt.close()


        # # XGBoost: OmniPath vs BioGRID
        # if 'xgboost' in preds_dict:
        #     preds_xgb = preds_dict['xgboost'][threshold]

        #     # Against BioGRID
        #     ref_biogrid = load_biogrid_reference(args.base_dir)
        #     preds_xgb_bio, ref_bio = filter_for_common_proteins(preds_xgb, ref_biogrid)
        #     merged_xgb_bio = pd.merge(
        #         preds_xgb_bio,
        #         ref_bio,
        #         on=['PredictiveFeature', 'TargetFeature'],
        #         how='left',
        #         indicator=True
        #     )
        #     y_true_xgb_bio = (merged_xgb_bio['_merge'] == 'both').astype(int).tolist()
        #     y_scores_xgb_bio = merged_xgb_bio[value_cols['xgboost']].tolist()
        #     fpr_xgb_bio, tpr_xgb_bio, _ = roc_curve(y_true_xgb_bio, y_scores_xgb_bio)
        #     roc_auc_xgb_bio = auc(fpr_xgb_bio, tpr_xgb_bio)

        #     # Against OmniPath
        #     ref_omnipath = load_omnipath_reference(args.base_dir)
        #     preds_xgb_omni, ref_omni = filter_for_common_proteins(preds_xgb, ref_omnipath)
        #     merged_xgb_omni = pd.merge(
        #         preds_xgb_omni,
        #         ref_omni,
        #         on=['PredictiveFeature', 'TargetFeature'],
        #         how='left',
        #         indicator=True
        #     )
        #     y_true_xgb_omni = (merged_xgb_omni['_merge'] == 'both').astype(int).tolist()
        #     y_scores_xgb_omni = merged_xgb_omni[value_cols['xgboost']].tolist()
        #     fpr_xgb_omni, tpr_xgb_omni, _ = roc_curve(y_true_xgb_omni, y_scores_xgb_omni)
        #     roc_auc_xgb_omni = auc(fpr_xgb_omni, tpr_xgb_omni)

        #     # Plot both ROC curves for XGBoost
        #     plt.figure(figsize=(7, 6))
        #     plt.plot(fpr_xgb_omni, tpr_xgb_omni, color='orange', lw=2, label=f'XGBoost OmniPath ROC (AUC = {roc_auc_xgb_omni:.2f})')
        #     plt.plot(fpr_xgb_bio, tpr_xgb_bio, color='green', lw=2, label=f'XGBoost BioGRID ROC (AUC = {roc_auc_xgb_bio:.2f})')
        #     plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        #     plt.xlim([0.0, 1.0])
        #     plt.ylim([0.0, 1.05])
        #     plt.xlabel('False Positive Rate', fontsize=14)
        #     plt.ylabel('True Positive Rate', fontsize=14)
        #     plt.title(f'XGBoost ROC: OmniPath vs BioGRID\n(Threshold {threshold}, SHAPValues)', fontsize=15)
        #     # plt.title(f'XGBoost ROC: OmniPath vs BioGRID\n(Threshold {threshold}, SHAP*R2)', fontsize=15)
        #     plt.legend(loc="lower right")
        #     plt.grid(True, linestyle='--', alpha=0.3)
        #     plt.tight_layout()
        #     plt.savefig(f"{args.base_dir}/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold{threshold}.png")
        #     # plt.savefig(f"{args.base_dir}/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold{threshold}_shap_r2.png")
        #     plt.close()


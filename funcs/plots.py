#!/bin/python

import pandas as pd
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
from pyvis.network import Network
from scipy.stats import norm
from sklearn.mixture import GaussianMixture


# ----------------- #

def scatterplot_phos_vals(dataset, save_file_as):
    plt.figure()
    # create a long-form dataframe
    dataset_melt = dataset.melt(id_vars='DatasetName', var_name='Phosphosite', value_name='Value')

    dataset_num = len(dataset)
    # create viridis colourmap 
    viridis = plt.get_cmap('viridis')
    colours = viridis(np.linspace(0, 0.85, dataset_num))
    colours = [mcolors.to_hex(c) for c in colours]

    # create scatter plot
    sns.stripplot(x='DatasetName', y='Value', hue='DatasetName', data=dataset_melt, jitter=True, 
                  size = 1.7, alpha = 0.9, 
                  palette=colours)

    # plot settings
    plt.xlabel('Dataset')
    plt.ylabel('Log2 abundance values')
    plt.title(f'Phosphoproteomics matrix log2-transformed data visualisation')
    plt.xticks(rotation=90)
    plt.grid(axis='y', linestyle='-')
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    
    plt.savefig(f'/data/home/bt24990/ExplainableAI/PlotsGraphsPNGs/{save_file_as}.png', dpi=300, bbox_inches='tight')
    print(f'Scatter plot of values saved successfully!')
    
    
# ----------------- #  

def lineplot_raw_phos_means(dataset, save_file_as):
    plt.figure()
    # plot mean value for each dataset
    means = dataset.set_index('DatasetName').mean(1)
    means.plot(kind='line')

    plt.xlabel('Dataset')
    plt.ylabel('Mean value')
    plt.title('Phosphoproteomics matrix mean dataset values')
    plt.xticks(rotation=90)

    # remove x-axis labels
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.savefig(f'/data/home/bt24990/ExplainableAI/PlotsGraphsPNGs/{save_file_as}.png', dpi=300, bbox_inches='tight')
    print('Line plot of mean values saved successfully!')
    
    

# ----------------- #  
    
def lineplot_norm_phos_means(dataset, save_file_as):
    plt.figure()
    # plot mean value for each dataset
    means = dataset.set_index('DatasetName').mean(1)
    means.plot(kind='line')

    plt.xlabel('Dataset')
    plt.ylabel('Mean value')
    plt.title('Phosphoproteomics matrix mean dataset values')
    plt.xticks(rotation=90)
    
    # Set the y scale between 0 and 1
    plt.ylim(0, 1)

    # remove x-axis labels
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.savefig(f'/data/home/bt24990/ExplainableAI/PlotsGraphsPNGs/{save_file_as}.png', dpi=300, bbox_inches='tight')
    print('Line plot of mean values saved successfully!')
    
    

# ----------------- #  

def plot_CNN_or_XGB_predicted_network(base_dir, model, df, strong_percentile, medium_percentile, 
                                      weak_percentile, selected_prots, save_as_filename):
    """Plots network of selected proteins.
    
    Inputs:
    df <dataframe>: 4 column df - PredictiveFeature, TargetFeature, SHAPValue,
        and LogSHAPValue columns
    strong_boundary <int>: value to set at strongest log SHAP value cutoff
    medium_boundary <int>: value to set at medium log SHAP value cutoff
    weak_boundary <int>: value to set at weakest log SHAP value cutoff 
    selected_prots <tuple>: Tuple of string proteins to select for
    save_as_filename <str>: Filename to save output network to
    """
    
    df['LogSHAPValue'] = pd.to_numeric(df['LogSHAPValue'], errors='coerce')
    abs_shap = df['LogSHAPValue'].abs().dropna()
    print('Absolute SHAP values:', abs_shap)
    strong_boundary = np.percentile(abs_shap, strong_percentile)
    medium_boundary = np.percentile(abs_shap, medium_percentile)
    weak_boundary = np.percentile(abs_shap, weak_percentile)

    strong_boundary = -strong_boundary
    medium_boundary = -medium_boundary
    weak_boundary = -weak_boundary

    print(f"Strong boundary: {strong_boundary}, Medium: {medium_boundary}, Weak: {weak_boundary}")
        
    nt = Network('1000px', '1000px', directed=True)
    for _, row in df.iterrows():
        PredictiveFeature = row['Feature']
        TargetFeature = row['TargetFeature']
        logSHAPvalue = row['LogSHAPValue']

        # Set default colors
        feature_color = '#a7b5e0'
        targetfeat_color = '#a7b5e0'

        # Update colors if features are in selected_prots
        if PredictiveFeature in selected_prots:
            feature_color = '#d16002'
        if TargetFeature in selected_prots:
            targetfeat_color = '#d16002'

        strong_color = '#311432'
        medium_color = '#a84296'
        weak_color = '#d8bfd8'

        node_fontsize = 25

        # for feedback loops
        if PredictiveFeature == TargetFeature:
            if logSHAPvalue >= strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=strong_color, value=2)
            elif logSHAPvalue >= medium_boundary and logSHAPvalue < strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=medium_color, value=2)
            elif logSHAPvalue >= weak_boundary and logSHAPvalue < medium_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=weak_color, value=2)
        # for all other edges
        else:
            if logSHAPvalue >= strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=strong_color, value=2)
            elif logSHAPvalue >= medium_boundary and logSHAPvalue < strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=medium_color, value=2)
            elif logSHAPvalue >= weak_boundary and logSHAPvalue < medium_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=weak_color, value=2)
            
    nt.save_graph(f'{base_dir}/08_results/{model}/predicted_networks/{save_as_filename}')
    print('Predicted network graph saved!')
    
    
# ----------------- #  

def plot_LR_predicted_network(base_dir, df, strong_percentile, medium_percentile, 
                              weak_percentile, selected_prots, save_as_filename):
    """Plots network of selected proteins.
    
    Inputs:
    df <dataframe>: 4 column df - PredictiveFeature, TargetFeature, SHAPValue,
        and LogSHAPValue columns
    strong_boundary <int>: value to set at strongest log SHAP value cutoff
    medium_boundary <int>: value to set at medium log SHAP value cutoff
    weak_boundary <int>: value to set at weakest log SHAP value cutoff 
    selected_prots <tuple>: Tuple of string proteins to select for
    save_as_filename <str>: Filename to save output network to
    """
    abs_coeff = df['MedianCoeff'].abs()
    strong_boundary = np.percentile(abs_coeff, 100 - strong_percentile)
    medium_boundary = np.percentile(abs_coeff, 100 - medium_percentile)
    weak_boundary = np.percentile(abs_coeff, 100 - weak_percentile)

    nt = Network('1000px', '1000px', directed=True)
    for row in df.iterrows():
        PredictiveFeature = row[1][0]
        TargetFeature = row[1][1]
        Coefficient = row[1][2]

        # Set default colors
        feature_color = '#a7b5e0'
        targetfeat_color = '#a7b5e0'

        # Update colors if features are in selected_prots
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
    
        # for feedback loops
        if PredictiveFeature == TargetFeature:
            if Coefficient <= -strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=strong_negative, value=2)
            elif Coefficient <= -medium_boundary and Coefficient > -strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=medium_negative, value=2)
            elif Coefficient <= -weak_boundary and Coefficient > -medium_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=weak_negative, value=2)
            elif Coefficient >= strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=strong_positive, value=2)
            elif Coefficient >= medium_boundary and Coefficient < strong_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=medium_positive, value=2)
            elif Coefficient >= weak_boundary and Coefficient < medium_boundary:
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, PredictiveFeature, color=weak_positive, value=2)
        # for all other edges
        else:
            if Coefficient <= -strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=strong_negative, value=2)
            elif Coefficient <= -medium_boundary and Coefficient > -strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=medium_negative, value=2)
            elif Coefficient <= -weak_boundary and Coefficient > -medium_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=weak_negative, value=2) 
            elif Coefficient >= strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=strong_positive, value=2)
            elif Coefficient >= medium_boundary and Coefficient < strong_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=medium_positive, value=2)
            elif Coefficient >= weak_boundary and Coefficient < medium_boundary:
                if TargetFeature not in nt.nodes:
                    nt.add_node(TargetFeature, color=targetfeat_color, font={'size': node_fontsize})
                if PredictiveFeature not in nt.nodes:
                    nt.add_node(PredictiveFeature, color=feature_color, font={'size': node_fontsize})
                nt.add_edge(PredictiveFeature, TargetFeature, color=weak_positive, value=2)
            
    nt.save_graph(f'{base_dir}/08_results/linear_regression/predicted_networks/{save_as_filename}')
    print('Predicted network graph saved!')
    
    
# ----------------- #  

def coefficients_histogram_for_specific_proteins(min_vals, proteins, network, strong_boundary, medium_boundary, weak_boundary):
    """Plots a histogram of coefficients against the number of predictive features for specified proteins."""
    
    df = pd.read_csv(f'/data/home/bt24990/ExplainableAI/LinearRegression/LR_Coefficients_Min{min_vals}Vals.csv', header=0) 

    # subset dataframe to only include rows with specified proteins
    dfs = []
    for i in proteins:
        target_feat_in_prots = df[df['TargetFeature'].str.contains(i)]
        dfs.append(target_feat_in_prots)
        predictive_feat_in_prots = df[df['Feature'].str.contains(i)]
        dfs.append(predictive_feat_in_prots)

    dfs = pd.concat(dfs)
    dfs = dfs[['TargetFeature', 'Feature', 'Coefficient']]
    sorted_df = dfs.sort_values(by=['Coefficient'], ascending=False)
    lr_coeff_vals = sorted_df['Coefficient']

    plt.figure()
    plt.hist(lr_coeff_vals, bins=100, log=True)
    plt.axvline(x=-strong_boundary, ls='--', color='#873c07')
    plt.axvline(x=-medium_boundary, ls='--', color='#ca5a0b')
    plt.axvline(x=-weak_boundary, ls='--', color='#f8bb8f')
    plt.axvline(x=strong_boundary, ls='--', color='#36441d')
    plt.axvline(x=medium_boundary, ls='--', color='#708d3e')
    plt.axvline(x=weak_boundary, ls='--', color='#d4e2bd')
    # plt.text(-6.8,10000,'strong\nnegative', fontsize=8, color='#873c07')
    # plt.text(-3.3, 10000,'medium\nnegative', fontsize=8, color='#ca5a0b')
    # plt.text(-2.3, 4500, 'weak\nnegative', fontsize=8, color='#f8bb8f')
    # plt.text(7.2, 10000,'strong\npositive', fontsize=8, color='#36441d')
    # plt.text(3.7, 4500,'medium\npositive', fontsize=8, color='#708d3e')
    # plt.text(2.7, 10000, 'weak\npositive', fontsize=8, color='#d4e2bd')
    plt.xlabel('Coefficient values')
    plt.ylabel('Frequency of features (log scale)')
    plt.title(f'Frequency distribution for LR \nCoefficient values ({network} proteins)')
    plt.savefig(f'/data/home/bt24990/ExplainableAI/PlotsGraphsPNGs/LR_CoefficientValues_{network}_Histogram_Min{min_vals}Vals.png', dpi=300, bbox_inches='tight')
    print(f"Histogram distribution of coefficients for Min{min_vals}Vals saved successfully!")
    

# ----------------- # 

def plot_logged_SHAP_values_for_each_threshold(model_type, network_name, 
                                                             log_shap_vals_50, log_shap_vals_100, 
                                                             log_shap_vals_150, log_shap_vals_200):
    """Plot 4 histograms (one figure) of logged SHAP values, one plot per threshold.
    
    This plot DOESN'T have distribution (Gaussian or normal) lines.
    """
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))
    data = [log_shap_vals_50, log_shap_vals_100, log_shap_vals_150, log_shap_vals_200]
    titles = ['Threshold: 50', 'Threshold: 100', 'Threshold: 150', 'Threshold: 200']
    for i, ax in enumerate(axs.flat):
        ax.hist(data[i], bins=100)
        ax.set_title(titles[i], fontsize=12)
        ax.set_xlabel('Logged SHAP values', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.grid(True, linestyle='-', alpha=0.3)
        # ax.grid(True, linestyle='-', alpha=0.3)
    # Set the overall figure title
    fig.suptitle(f'Comparison of {model_type} logged SHAP values across\ndifferent thresholds (models optimised for MSE)')
    fig.subplots_adjust(hspace=0.35)
    plt.savefig(f'/data/home/bt24990/ExplainableAI/07_results/{model_type}/figures/{model_type}_loggedSHAPvalues_{network_name}_histogram_all_thresholds.png', dpi=300, bbox_inches='tight')
    
    
# ----------------- #  
    
def plot_logged_SHAP_values_distribution_for_one_threshold(min_vals, model_type, network_name, log_shap_vals):
    """Plot a histogram of logged SHAP values for a specific threshold."""
        
    gmm = GaussianMixture(n_components=2)
    log_shap_vals_reshaped = log_shap_vals.reshape(-1, 1)
    gmm.fit(log_shap_vals_reshaped)
    mu1, mu2 = gmm.means_.flatten()
    std1, std2 = np.sqrt(gmm.covariances_).flatten()
    weight1, weight2 = gmm.weights_
    
    plt.figure()
    count, bins, ignored = plt.hist(log_shap_vals, bins=100, density=False)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p1 = norm.pdf(x, mu1, std1)
    p2 = norm.pdf(x, mu2, std2)
    p1_scaled = p1 * weight1 * len(log_shap_vals) * (bins[1] - bins[0])
    p2_scaled = p2 * weight2 * len(log_shap_vals) * (bins[1] - bins[0])
    p_combined = p1_scaled + p2_scaled

    # combine the two normal distributions
    plt.plot(x, p_combined, 'k', linewidth=2)
    plt.xlabel('Logged SHAP values')
    plt.ylabel('Frequency')
    plt.grid(True, linestyle='-', alpha=0.3)
    plt.title(f'Distribution for {model_type} logged SHAP values\n({network_name} proteins) - threshold: {min_vals}, optimised for MSE')
    plt.savefig(f'/data/home/bt24990/ExplainableAI/07_results/{model_type}/figures/{model_type}_loggedSHAPvalues_{network_name}_histogram_min{min_vals}vals.png', dpi=300, bbox_inches='tight')
    
    print(f'Plot of logged CNN SHAP values (threshold: {min_vals}) saved successfully!')
    
    
# ----------------- #

def plot_logged_SHAP_values_distribution_for_all_thresholds(model_type, network_name, 
                                                             log_shap_vals_50, log_shap_vals_100, 
                                                             log_shap_vals_150, log_shap_vals_200):
    """Plot 4 individual histograms of logged SHAP values, one plot per threshold.
    
    distributions (list): List of distribution types for each dataset. Options: 'normal', 'gaussian', or None.
    """
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
    data = [log_shap_vals_50, log_shap_vals_100, log_shap_vals_150, log_shap_vals_200]
    titles = ['Threshold: 50', 'Threshold: 100', 'Threshold: 150', 'Threshold: 200']
    
    for i, ax in enumerate(axs.flat):
        # Plot the histogram
        count, bins, _ = ax.hist(data[i], bins=100, density=False)

        gmm = GaussianMixture(n_components=2)
        reshaped_data = np.array(data[i]).reshape(-1, 1)
        gmm.fit(reshaped_data)
        mu1, mu2 = gmm.means_.flatten()
        std1, std2 = np.sqrt(gmm.covariances_).flatten()
        weight1, weight2 = gmm.weights_
        
        xmin, xmax = ax.get_xlim()
        x = np.linspace(xmin, xmax, 100)
        bin_width = bins[1] - bins[0]
        p1 = norm.pdf(x, mu1, std1) * weight1
        p2 = norm.pdf(x, mu2, std2) * weight2
        p_combined = (p1 + p2) * bin_width * len(data[i])
        ax.plot(x, p_combined, 'k', linewidth=2, label=f'GMM: μ1={mu1:.2f}, μ2={mu2:.2f}')
        
        # Add labels and titles
        ax.set_title(titles[i], fontsize=12)
        ax.set_xlabel('Logged SHAP values', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.grid(True, linestyle='-', alpha=0.3)
    
    # Set the overall figure title
    fig.suptitle(f'Comparison of {model_type} logged SHAP values for {network_name}\nnetwork across different thresholds (models optimised for MSE)', fontsize=16)
    fig.subplots_adjust(hspace=0.4, wspace=0.3)
    plt.savefig(f'/data/home/bt24990/ExplainableAI/07_results/{model_type}/figures/{model_type}_loggedSHAPvalues_{network_name}_histogram_all_thresholds.png', dpi=300, bbox_inches='tight')
    plt.show()
    
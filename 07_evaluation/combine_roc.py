import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

def combine_pngs(plot_paths, out_path, nrows=2, ncols=4, figsize=(32, 16)):
    fig, axs = plt.subplots(nrows, ncols, figsize=figsize)
    for i, ax in enumerate(axs.flat):
        if i < len(plot_paths) and os.path.exists(plot_paths[i]):
            img = mpimg.imread(plot_paths[i])
            ax.imshow(img)
        ax.axis('off')
    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.01, hspace=0.01)
    fig.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0)
    plt.show()

# Linear Regression plots
lr_plot_paths = [
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold50.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold100.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold150.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold200.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold50_coeff_r2.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold100_coeff_r2.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold150_coeff_r2.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_threshold200_coeff_r2.png",
]
combine_pngs(
    lr_plot_paths,
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/LR_omnipath_biogrid_roc_curve_all_thresholds.png"
)

# XGBoost plots
xgb_plot_paths = [
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold50.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold100.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold150.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold200.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold50_shap_r2.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold100_shap_r2.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold150_shap_r2.png",
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGBoost_omnipath_biogrid_roc_curve_threshold200_shap_r2.png",
]
combine_pngs(
    xgb_plot_paths,
    "/data/home/bt24990/ExplainableAI/08_results/roc_plots/XGB_omnipath_biogrid_roc_curve_all_thresholds.png"
)
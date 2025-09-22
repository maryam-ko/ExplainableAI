import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os

def combine_pngs(plot_paths, out_path, nrows=1, ncols=3, figsize=(20, 8)):
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

# omnipath plots
omnipath_plot_paths = [
    "/data/home/bt24990/ExplainableAI/08_results/omnipath_results/omnipath_aucs_all_interactions.png",
    "/data/home/bt24990/ExplainableAI/08_results/omnipath_results/omnipath_aucs_conserved_by_model.png",
    "/data/home/bt24990/ExplainableAI/08_results/omnipath_results/omnipath_aucs_conserved_by_threshold.png",
]
combine_pngs(
    omnipath_plot_paths,
    "/data/home/bt24990/ExplainableAI/08_results/omnipath_results/omnipath_combined.png"
)
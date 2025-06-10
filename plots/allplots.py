import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Define the image paths and titles in the desired order
plot_paths = [
    "/data/home/bt24990/ExplainableAI/plots/raw_scatter.png",
    "/data/home/bt24990/ExplainableAI/plots/normalise_scatter.png",
    "/data/home/bt24990/ExplainableAI/plots/arrays_scatter.png",
    "/data/home/bt24990/ExplainableAI/plots/Zscore_scatter.png",
    "/data/home/bt24990/ExplainableAI/plots/raw_mean_plot.png",
    "/data/home/bt24990/ExplainableAI/plots/normalised_mean_plot.png",
    "/data/home/bt24990/ExplainableAI/plots/arrays_mean_plot.png",
    "/data/home/bt24990/ExplainableAI/plots/Zscore_mean_plot.png"
]

# Create a 2x3 panel layout
fig, axs = plt.subplots(2, 4, figsize=(30, 15))

for i, ax in enumerate(axs.flat):
    img = mpimg.imread(plot_paths[i])
    ax.imshow(img, aspect='equal', interpolation='none')
    ax.axis('off')

plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.01, hspace=0.01)

# Tidy up layout
fig.tight_layout()

# Show or save
plt.show()
plt.savefig("matrix_plots.png", dpi=300, bbox_inches='tight', pad_inches=0)

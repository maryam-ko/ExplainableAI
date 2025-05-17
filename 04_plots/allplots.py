import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Define the image paths and titles in the desired order
plot_paths = [
    "/data/home/bt24990/maryam-ko-QMUL-MSc-Project/04_plots/raw_scatter_plot.png",
    "/data/home/bt24990/maryam-ko-QMUL-MSc-Project/04_plots/normalised_scatter.png",
    "/data/home/bt24990/maryam-ko-QMUL-MSc-Project/04_plots/raw_mean_plot.png",
    "/data/home/bt24990/maryam-ko-QMUL-MSc-Project/04_plots/normalised_mean.png"
]

# Create a 2x2 panel layout
fig, axs = plt.subplots(2, 2, figsize=(12, 10))

# Loop and assign images
for i, ax in enumerate(axs.flat):
    img = mpimg.imread(plot_paths[i])
    ax.imshow(img)
    ax.axis('off') 

# Tidy up layout
plt.tight_layout()

# Show or save
plt.show()
plt.savefig("comparison_panels.png", dpi=300)

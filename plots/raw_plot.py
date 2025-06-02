import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load your dataset
df = pd.read_csv("/data/home/bt24990/maryam-ko-QMUL-MSc-Project/02_raw_matrix/MatrixCSVs/RawMatrix_NoOutliers.csv")

# Melt the dataframe for plotting (long format)
df_melted = df.melt(id_vars='DatasetName', var_name='PhosphoSite', value_name='Quantification Value')

dataset_order = df['DatasetName'].drop_duplicates().reset_index(drop=True)
dataset_index = {name: i for i, name in enumerate(dataset_order)}
df_melted['ColorIndex'] = df_melted['DatasetName'].map(dataset_index)

norm = plt.Normalize(df_melted['ColorIndex'].min(), df_melted['ColorIndex'].max())
colors = plt.cm.magma(norm(df_melted['ColorIndex']))  # You can change 'plasma' to other colormaps like 'viridis', 'inferno', etc.

# Plot 1: Raw Scatter Plot of Quantification Values
plt.figure(figsize=(7, 5))
plt.scatter(
    x=df_melted['DatasetName'],
    y=df_melted['Quantification Value'],
    c=colors,
    alpha=0.5,
    s=5
)

plt.xticks([], [])  # Hides x-axis ticks and labels
plt.xlabel("DatasetName")
plt.ylabel("Quantification Value")
plt.title("Raw log2 phosphorylation values\n"
          "(before normalisation)")
plt.tight_layout(pad=2.0)  # Increase padding between title and plot
plt.savefig("raw_scatter.png", dpi=300)
plt.close()

# Plot 2: Mean Quantification Value per Dataset
df_means = df.drop(columns='DatasetName').mean(axis=1)
plt.figure(figsize=(7, 5))
plt.plot(df['DatasetName'], df_means, color='indigo')
plt.xticks([])  
plt.xlabel("DatasetName")
plt.ylabel("Mean Quantification Value")
plt.tight_layout()
plt.savefig("raw_mean_plot.png", dpi=300)
plt.close()

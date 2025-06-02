import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load your dataset
df = pd.read_csv("/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NormalisedMatrix.csv")
# Check duplicates
print("Duplicates:", df.columns[df.columns.duplicated()])

# Remove duplicates
df = df.loc[:, ~df.columns.duplicated()]

# Rename index column if present (sometimes 'index' is a column)
if 'index' in df.columns:
    df = df.rename(columns={'index': 'DatasetName'})

# Check 'DatasetName' exists exactly once
assert df.columns.tolist().count('DatasetName') == 1, "Multiple DatasetName columns!"

# Now melt
df_melted = df.melt(id_vars='DatasetName', var_name='PhosphoSite', value_name='Quantification Value')

# Map DatasetName to numeric index
dataset_order = df['DatasetName'].drop_duplicates().reset_index(drop=True)
dataset_index = {name: i for i, name in enumerate(dataset_order)}
df_melted['DatasetIndex'] = df_melted['DatasetName'].map(dataset_index)

# Normalize for color mapping
norm = plt.Normalize(df_melted['DatasetIndex'].min(), df_melted['DatasetIndex'].max())
colors = plt.cm.magma(norm(df_melted['DatasetIndex']))

# Plot: Scatter with cubehelix gradient
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
plt.title("Min-max scaled log2 phosphorylation values\n"
          "(normalised per dataset)")
plt.tight_layout(pad=2.0)  # Increase padding between title and plot
plt.savefig("normalise_scatter.png", dpi=300)
plt.close()

# Plot 2: Mean Quantification Value per Dataset
df_means = df.drop(columns='DatasetName').mean(axis=1)

plt.figure(figsize=(7, 5))
plt.plot(df['DatasetName'], df_means, color='indigo')

# Hide x-axis labels and ticks
plt.xticks([])  
plt.xlabel("DatasetName")
plt.ylabel("Mean Quantification Value")
plt.tight_layout()

# Save the figure
plt.savefig("normalised_mean_plot.png", dpi=300)
plt.close()

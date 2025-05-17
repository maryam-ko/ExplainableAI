import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your dataset
df = pd.read_csv("/data/home/bt24990/maryam-ko-QMUL-MSc-Project/03_normalise_data/MatrixCSVs/NormalisedMatrix.csv")

# Melt the dataframe for plotting (long format)
df_melted = df.melt(id_vars='DatasetName', var_name='PhosphoSite', value_name='Quantification Value')

# Plot 1: Raw Scatter Plot of Quantification Values
plt.figure(figsize=(10, 6))
ax = sns.stripplot(
    data=df_melted,
    x='DatasetName',
    y='Quantification Value',
    color='orange',
    alpha=0.5,
    jitter=True,
    size=3
)

# Hide the x-axis ticks and labels
ax.set_xticks([])  # Hides the x-axis ticks
ax.set_xticklabels([])  # Hides the x-axis labels

plt.title("Min-max scaled log2 phosphorylation values (normarlised per dataset)")
plt.xlabel("DatasetName")
plt.ylabel("Quantification Value")
plt.tight_layout()

# Save the figure
plt.savefig("normalised_scatter_clean.png", dpi=300)
plt.close()

# Plot 2: Mean Quantification Value per Dataset
df_means = df.drop(columns='DatasetName').mean(axis=1)

plt.figure(figsize=(10, 6))
plt.plot(df['DatasetName'], df_means, color='dodgerblue')

# Hide x-axis labels and ticks
plt.xticks([])  
plt.xlabel("DatasetName")
plt.ylabel("Mean Quantification Value")
plt.tight_layout()

# Save the figure
plt.savefig("normalised_mean_plot.png", dpi=300)
plt.close()

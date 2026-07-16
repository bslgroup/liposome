import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set(style="white") 
#sns.set(style="whitegrid")

directory = "./"

file_names = ["0_sph_leaflet_data.txt", "10_sph_leaflet_data.txt", "20_sph_leaflet_data.txt", "30_sph_leaflet_data.txt", "40_sph_leaflet_data.txt"]

column_names = ['OuterRadius', 'InnerRadius', 'OuterStd Dev', 'InnerStdDev']
x_labels = ["Radius (Å)", "Radius (Å)", "Outer Leaflet Smoothness (Distance Stdev)", "Inner Leaflet Smoothness (Distance Stdev)"]
colors = ['blue', 'orange', 'green', 'red', 'purple']
labels = [name.split('_')[0] for name in file_names]

plt.figure(figsize=(10, 8))
for j, file_name in enumerate(file_names):
    file_path = os.path.join(directory, file_name)
    df = pd.read_csv(file_path, sep="\t")
    sns.kdeplot(df['OuterRadius'], label=f"{labels[j]} Outer", color=colors[j], fill=True, alpha=0.5)
    sns.kdeplot(df['InnerRadius'], label=f"{labels[j]} Inner", color=colors[j], linestyle="--", fill=True, alpha=0.3)
plt.title("Probability Density of Radius")
plt.xlabel("Radius (Å)")
plt.ylabel("Probability Density")
plt.legend()
plt.savefig("RadiusDistribution.pdf", dpi=300)
plt.clf()

def plot_density(column, x_label, file_name):
    plt.figure(figsize=(10, 8))
    for j, fname in enumerate(file_names):
        path = os.path.join(directory, fname)
        df = pd.read_csv(path, sep="\t")
        sns.kdeplot(df[column], label=labels[j], color=colors[j], fill=True, alpha=0.5)
    plt.title(f"Probability Density of {x_label}")
    plt.xlabel(x_label)
    plt.ylabel("Probability Density")
    plt.legend()
    plt.savefig(file_name, dpi=300)
    plt.clf()

plot_density('OuterStd Dev', x_labels[2], "OuterLeafletSmoothness.pdf")
plot_density('InnerStdDev', x_labels[3], "InnerLeafletSmoothness.pdf")


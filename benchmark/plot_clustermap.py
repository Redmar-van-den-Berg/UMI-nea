import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import dendrogram, linkage

infile = "umi_count_7_clonotypes.csv"
df = pd.read_csv(infile)
df['feature'] = (
    df['platform'] + '_' + 
    df['condition'].astype(str) + '_' + 
    df['replicate'] + '_' +
    df['clonotype'].astype(str)
)
tools_ord = ["Calib","UMIc-seq","UMI-tools","UMI-nea", "MiXCR"]
df['tools'] = pd.Categorical(df['tools'],categories=tools_ord,ordered=True)
df = df.sort_values(by=["tools"])
tool_matrix = df.pivot(index='tools', columns='feature', values='number of UMI')
corr_matrix = tool_matrix.T.corr()
mask = np.ones_like(corr_matrix, dtype=bool)
for i in range(0,5):
    for j in range(0,5):
        if i<j:
            mask[i,j] = False
        if i==j:
            mask[i,j] = False

plt.figure(figsize=(10, 8))
g = sns.clustermap(
    corr_matrix.iloc[[1,3,4,0,2], :].iloc[:, [1,3,4,0,2]],
    cmap='coolwarm',
    annot=True,
    fmt=".3f",
    row_cluster=False, col_cluster=False,
    mask = mask,
    annot_kws={"size": 15},
    vmin=0.84,
    vmax=1,
    linewidths=1,
)
divider = make_axes_locatable(g.ax_heatmap)
ax_row_dendro = divider.append_axes("top", size="20%", pad=0.1)
Z_row = linkage(corr_matrix, method='ward')
dendrogram(
    Z_row,
    orientation='top',
    ax=ax_row_dendro,
    color_threshold=0
)
ax_row_dendro.axis('off')
g.ax_heatmap.set_xlabel("Tools", fontsize=15)
g.ax_heatmap.set_ylabel("Tools", fontsize=15)
g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(), 
    fontsize=14)
g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_yticklabels(), 
    fontsize=14)
cb = g.cax
cb.tick_params(labelsize=14)
plt.savefig("clonotype_umi_count_linkage.png", bbox_inches="tight", dpi=350)
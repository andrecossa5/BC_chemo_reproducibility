"""
Script to analyze pyscenic results
"""

from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose


import os
import json
import base64
import zlib
import pandas as pd 
import seaborn as sns
import numpy as np
import scanpy as sc
import loompy as lp
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss


from IPython.display import HTML, display


#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility"
path_adata = '/Users/ieo7295/Desktop/tests/cell/data/default'
path_data= os.path.join(path_main,"data", "CTCs","grn")
path_results=os.path.join(path_main,"results", "CTCs", "gene_reg_net")


#Data
adata = sc.read_h5ad(os.path.join(path_adata, 'clustered.h5ad'))
lf = lp.connect(os.path.join(path_data, "pyscenic_res", "ctc_aucell_pruned.loom") , mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons= lf.ra.Regulons
lf.close()

#Extract metadata

meta = adata.obs[['GBC', 'sample', 'dataset', 'origin', 'nUMIs',
                   'mito_perc', 'detected_genes', 'cell_complexity',
                     'n_genes', 'seq_run','100_NN_1.5', 'leiden', 'final_cell_state']]
lf = lp.connect(os.path.join(path_data, "pyscenic_res", "scenic_visualize.loom"), mode='r+', validate=False )
sc_mtx= pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

lf.ra.keys()
lf.ca.keys()
lf.attrs.keys()

#Extract thresholds
metadata_dict = json.loads(lf.attrs["MetaData"]) 
metadata_dict['regulonThresholds']
filtered_meta_dict= {item['regulon']: item['defaultThresholdValue'] for item in metadata_dict['regulonThresholds']}
thresholds = {k.replace('_', ''): v for k, v in filtered_meta_dict.items()}
threshold_series = pd.Series(thresholds).reindex(sc_mtx.columns)
threshold_bin= {key.rstrip("(+)"): value for key, value in thresholds.items()}
sc_mtx.rename(columns=lambda col: col.rstrip("(+)"), inplace=True)

#Binarization
fig, (ax1, ax2, ax3,) = plt.subplots(1, 3, figsize=(8, 4), dpi=100)

plot_binarization(sc_mtx, 'ZNF384', threshold_bin['ZNF384'], ax=ax1)
plot_binarization(sc_mtx, 'ARID5B', threshold_bin['ARID5B'], ax=ax2)
plot_binarization(sc_mtx, 'ATF3', threshold_bin['ATF3'], ax=ax3)

plt.tight_layout()
plt.show()


#Heatmap with binarized regulon activity
def palplot(pal, names, colors=None, size=1):
    n = len(pal)
    f, ax = plt.subplots(1, 1, figsize=(n * size, size))
    ax.imshow(np.arange(n).reshape(1, n),
              cmap=mpl.colors.ListedColormap(list(pal)),
              interpolation="nearest", aspect="auto")
    ax.set_xticks(np.arange(n) - .5)
    ax.set_yticks([-.5, .5])
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    colors = n * ['k'] if colors is None else colors
    for idx, (name, color) in enumerate(zip(names, colors)):
        ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
    return f

N_COLORS = len(adata.obs['final_cell_state'].dtype.categories)
COLORS = [color['color'] for color in mpl.rcParams['axes.prop_cycle']]

cell_type_color_lut = dict(zip(adata.obs['final_cell_state'].dtype.categories, COLORS))
cell_id2cell_type_lut = meta.rename_axis('cell_id')['final_cell_state'].to_dict()
bw_palette = sns.xkcd_palette(["white", "black"])

sns.set()
sns.set_style("whitegrid")
fig = palplot(bw_palette, ['OFF', 'ON'], ['k', 'w'])
plt.show()

sns.set()
sns.set(font_scale=1.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
g = sns.clustermap(sc_mtx.T, 
               col_colors=sc_mtx.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
               cmap=bw_palette, figsize=(20,20))
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cell Type'])
g.cax.set_visible(False)
plt.show()
#g.fig.savefig(os.path.join(FIGURES_FOLDERNAME, 'clustermap - GSE115978.png'), format='png')

col_colors=sc_mtx.index.map(cell_id2cell_type_lut).map(cell_type_color_lut)
"""
Script to perform visualization of pseudobulk DE.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import scanpy as sc
from anndata import AnnData
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'pseudobulk')

# Read clustered
clustered = sc.read(os.path.join(path_data, 'clustered.h5ad'))
clustered.obs['GBC_sample'] = clustered.obs['GBC'].astype(str) + clustered.obs['sample'].astype(str)
GBC_size = clustered.obs.groupby('GBC_sample').size().loc[lambda x: x>=10]

# Read aggregated clone table
df = pd.read_csv(os.path.join(path_data, 'agg_pseudobulk.csv'), index_col=0)
meta_columns = df.columns[:8].to_list()

# Create AnnData object
adata = AnnData(
    X=df.loc[:, ~df.columns.isin(meta_columns)].values, 
    obs=df[meta_columns].join(clustered.obs.groupby('GBC_sample').size().to_frame('size')), 
    var=pd.DataFrame(index=df.columns[~df.columns.isin(meta_columns)])
)

# PP, PCA, UMAP
sc.pp.normalize_total(adata, target_sum=10**6)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', inplace=True, subset=True)
sc.pp.scale(adata)
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=10)
sc.tl.umap(adata)
df_embs = adata.obs.join(
    pd.DataFrame(
        adata.obsm['X_pca'], index=adata.obs_names,
        columns=[ f'PC{i+1}' for i in range(adata.obsm['X_pca'].shape[1])],
    )
    .join(
    pd.DataFrame(
        adata.obsm['X_umap'], index=adata.obs_names,
        columns=[ f'UMAP{i+1}' for i in range(adata.obsm['X_umap'].shape[1])],
    )
))


# Viz 

fig, ax = plt.subplots(figsize=(4,4))
colors = create_palette(df_embs, 'seq_run', 'tab10')
scatter(
    df_embs, x='UMAP1', y='UMAP2', c=colors, 
    by='seq_run', a=.5, s='size', scale_x=.7, ax=ax
)
format_ax(ax, title='Lentiviral clones')
add_legend(label='Seq run', colors=colors, ax=ax, loc='lower left', 
           bbox_to_anchor=(0,0), label_size=8, artists_size=6, ticks_size=6)
ax.axis('off')
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(5,5))
colors = create_palette(df_embs, 'condition', 'tab10')
scatter(
    df_embs, x='UMAP1', y='UMAP2', c=colors, 
    by='condition', a=.5, s='size', scale_x=.7, ax=ax
)
format_ax(ax, title='Lentiviral clones')
add_legend(label='Condition', colors=colors, ax=ax, loc='lower left', 
           bbox_to_anchor=(0,0), label_size=10, artists_size=8, ticks_size=8)
ax.axis('off')
fig.tight_layout()
plt.show()





# Visualization




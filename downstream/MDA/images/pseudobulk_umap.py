"""
UMAP of pseudobulk clones.
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
clustered.obs['GBC_sample'] = clustered.obs['GBC'].astype(str) + \
                              clustered.obs['sample'].astype(str)
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
seed = 0
sc.pp.normalize_total(adata, target_sum=10**6)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='cell_ranger', inplace=True, subset=True)
sc.pp.scale(adata)
sc.pp.pca(adata, n_comps=30, random_state=seed)
sc.pp.neighbors(adata, n_neighbors=10, random_state=seed)
sc.tl.umap(adata, random_state=seed)
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

# Read colors
clustered.obs['mock'] = clustered.obs['GBC'].astype(str) + clustered.obs['sample'].astype(str)
df_embs = df_embs.drop(columns=['condition'])
df_embs['condition'] = df_embs.index.map(
    lambda x: clustered.obs.loc[clustered.obs['mock']==x, 'condition'].unique()[0]
)
df_embs['condition'] = pd.Categorical(
    df_embs['condition'], 
    categories=clustered.obs['condition'].cat.categories
)
condition_colors = create_palette(df_embs, 'condition', ten_godisnot)

# Save
df_embs.to_csv(os.path.join(path_data, 'embs_pseudo.csv'))


##


# Viz 
fig, ax = plt.subplots(figsize=(6.5,5))

scatter(
    df_embs, x='UMAP1', y='UMAP2', c=condition_colors, 
    by='condition', a=.5, s='size', scale_x=.7, ax=ax
)
format_ax(ax, title='Lentiviral clones')
add_legend(label='Condition', colors=condition_colors, ax=ax, loc='upper left', 
           bbox_to_anchor=(1,1), label_size=10, artists_size=8, ticks_size=8)
ax.axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'umap_pseudobulk.png'), dpi=400)


##




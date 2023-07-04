"""
Make a custom embedding. 
"""

import os
import re
import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from anndata import AnnData
from Cellula.preprocessing._pp_recipes import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._embeddings import *
from Cellula.preprocessing._neighbors import *
from Cellula.preprocessing._integration import *
from Cellula.clustering._clustering import *
from Cellula.plotting._colors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'embeddings')


##


############################### PP and visualization, no_cc, Harmony

# Read
adata = sc.read(os.path.join(path_data, 'subsampled.h5ad'))
adata.raw = adata.copy()
adata.obs = adata.obs.iloc[:,:10]

# Summary sample QC
(
    adata.obs
    .loc[:, ['nUMIs', 'detected_genes', 'mito_perc', 'sample']]
    .groupby('sample')
    .agg('median')
)

# PP
name = 'regress_cc_5000'

_, reduced = regress_cc_pp_recipe(adata, n_HVGs=5000)
pca(reduced, layer='regressed', auto=True, GSEA=True, 
    return_adata=True, path_viz=path_results)

# compute_scVI(
#     reduced, 
#     n_latent=30, 
#     max_epochs=600,
#     layer='raw', 
#     categorical='study', 
#     continuous=['mito_perc', 'nUMIs', 'cycling', 'ribo_genes']
# )
compute_kNN(reduced, layer='regressed', int_method='original', k=15)
leiden(reduced, obsp_key='regressed|original|X_pca|NN_conn', res=.5)
leiden(reduced, obsp_key='regressed|original|X_pca|NN_conn', 
    obs_key='leiden_high_resolution', res=3)
df = embeddings(reduced, layer='regressed', rep='original', 
    paga_groups='leiden_high_resolution', umap_only=False)

# Save full embs and complete adata with pp slots
# df.to_csv(os.path.join(path_data, 'embeddings.full'))
# _.obsm['X_pca'] = reduced.obsm['regressed|original|X_pca']
# _.obsm['X_emb'] = df.loc[:, ['UMAP1', 'UMAP2']].values
# _.obsm['X_scVI'] = reduced.obsm['original|scVI|X_corrected']
# _.obs['leiden'] = df['leiden']
# _.obsp = reduced.obsp
# sanitize_neighbors(_, obsm_key='X_pca', old_neighbors_key='raw|scVI|X_corrected')
# _.write(os.path.join(path_data, 'subsampled_with_embs.h5ad'))

# Read
# df = pd.read_csv(os.path.join(path_data, 'embeddings.full'), index_col=0)
# adata = sc.read(os.path.join(path_data, 'subsampled_with_embs.h5ad'))
df['leiden'] = pd.Categorical(df['leiden'])
df['leiden_high_resolution'] = pd.Categorical(df['leiden_high_resolution'])

##

# Top clones
top_clones = df.groupby('GBC').size().loc[lambda x: x.apply(lambda x: x>10)].index
df['top_clones'] = df['GBC'].map(lambda x: x if x in top_clones else 'others')

# Plot cat
fig, axs = plt.subplots(2,3,figsize=(10,7))
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='leiden', ax=axs[0,0], s=4)
                #axes_kwargs={'legend':False})
draw_embeddings(df, 'FA1', 'FA2', cat='top_clones', ax=axs[0,1], s=4,
                axes_kwargs={'legend':False})
draw_embeddings(df, 'FA1', 'FA2', cat='sample', ax=axs[0,2], s=4, 
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='seq_run', ax=axs[1,0], s=4)
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='origin', ax=axs[1,1], s=4)
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='leiden_high_resolution', ax=axs[1,2], s=4,
                axes_kwargs={'legend':False})
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, f'{name}_cat_embeddings.png'
    ), 
    dpi=500
)

# Plot cont
fig, axs = plt.subplots(1,4,figsize=(15,3.5))
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='nUMIs', ax=axs[0], s=2,
    cbar_kwargs={'vmin':800, 'vmax':7500, 'layout':( (1.05,.25,.03,.5), 'right' )}
)
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='cycling', ax=axs[1], s=2,
    cbar_kwargs={'vmin':-.1, 'vmax':.1, 'layout':( (1.05,.25,.03,.5), 'right' )}
)
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='ribo_genes', ax=axs[2], s=2,
    cbar_kwargs={'vmin':0, 'vmax':2.5, 'layout':( (1.05,.25,.03,.5), 'right' )}
)
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='mito_perc', ax=axs[3], s=2,
    cbar_kwargs={'vmin':0, 'vmax':.12, 'layout':( (1.05,.25,.03,.5), 'right' )}
)
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, f'{name}_cont_embeddings.png'
    ), 
    dpi=500
)

# Biplots
df_ = (
    df.loc[:, df.columns[df.columns.str.contains('Diff')]].iloc[:,:5]
    .join(reduced.obs.loc[:, ['origin', 'seq_run', 'mito_perc', 'nUMIs', 'leiden']])
)
colors = None

# Biplot Diff comp
diff_comp = df_.columns[:-5]
from itertools import combinations
combos = [ (x, y) for x,y in combinations(diff_comp, 2) if x != y ]
make_folder(path_results, f'{name}_diff_comp')


for combo in combos:
    fig, axs = plt.subplots(1,5,figsize=(15.5, 3))
    x, y = combo
    draw_embeddings(df_, x, y, cat='origin', ax=axs[0], s=4)
    draw_embeddings(df_, x, y, cat='seq_run', ax=axs[1], s=4)
    draw_embeddings(df_, x, y, cont='mito_perc', ax=axs[2], s=4,
        cbar_kwargs={'vmin':0, 'vmax':.12, 'layout':( (1.05,.25,.03,.5), 'right' )})
    draw_embeddings(df_, x, y, cont='nUMIs', ax=axs[3], s=4,
        cbar_kwargs={'vmin':500, 'vmax':15000, 'layout': 'v2'})
    draw_embeddings(df_, x, y, cat='leiden', ax=axs[4], s=4)
    fig.tight_layout()
    fig.savefig(
    os.path.join(
        path_results, f'{name}_diff_comp', f'{x}_{y}.png'
    ), 
    dpi=500
)
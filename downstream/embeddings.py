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


#


# Utils
def remove_unwanted(a):
    
    mt = set(a.var_names[a.var_names.str.startswith('MT-')])
    ribo_1 = set(a.var_names[a.var_names.str.startswith('RPL')])
    ribo_2 = set(a.var_names[a.var_names.str.startswith('RPS')])
    anti_1 = set(a.var_names[a.var_names.str.contains('-AS')])
    anti_2 = set(a.var_names[a.var_names.str.startswith('AC0')])
    to_exclude = mt | ribo_1 | ribo_2 | anti_1 | anti_2
    a = a[:, ~a.var_names.isin(to_exclude)].copy()
    
    return a


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'embeddings')

##

# Read 
adata_original = sc.read(os.path.join(path_data, 'subsampled.h5ad'))

##

######################################### Default, no Harmony
a = AnnData(
    adata_original.layers['raw'], 
    obs=adata_original.obs.iloc[:,:11], 
    var=adata_original.var.iloc[:, [0,1]]
)
a.raw = a.copy()

# Remove MT-, RB-, noncoding genes.

# MOLTO BELLO 5000 hvgs, auto pcs, 30 k, 1 res, FA_diff.
_, reduced = standard_pp_recipe(a, n_HVGs=5000)
reduced = remove_unwanted(reduced)
make_folder(path_results, 'standard_scaled')
d_pca = pca(reduced, layer='scaled', auto=True, GSEA=True, path_viz=os.path.join(path_results, 'standard_scaled'))
reduced.obsm['scaled|original|X_pca'] = d_pca['X_pca'] 
compute_kNN(reduced, layer='scaled', k=30)
reduced.obs['leiden_high_res'] = leiden_clustering(
    reduced.obsp['scaled|original|X_pca|NN_conn'], res=1
)
reduced.obs['leiden_high_res'] = pd.Categorical(reduced.obs['leiden_high_res']) # High res for embeddings
reduced.obs['leiden'] = leiden_clustering(
    reduced.obsp['scaled|original|X_pca|NN_conn'], res=.5                       # Low res for clustering itself
)
reduced.obs['leiden'] = pd.Categorical(reduced.obs['leiden'])

# Embs
df = embeddings(reduced, paga_groups='leiden_high_res', umap_only=False)

# Cross
pd.crosstab(df['leiden'], df['seq_run'])
pd.crosstab(df['origin'], df['seq_run'])
pd.crosstab(df['condition'], df['seq_run'])

# fig, ax=plt.subplots()
# plot_heatmap(pd.crosstab(df['condition'], df['seq_run']), ax=ax)
# plt.show()
# N.B. Good condition-seq-run corresponance

# Top clones
top_clones = df.groupby('GBC').size().loc[lambda x: x.apply(lambda x: x>50)].index
df['top_clones'] = df['GBC'].map(lambda x: x if x in top_clones else 'others')

# Plot
fig, axs = plt.subplots(2,3,figsize=(10,7))
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='leiden', ax=axs[0,0], s=2)
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='top_clones', ax=axs[0,1], s=2,
                axes_kwargs={'legend':False})
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='sample', ax=axs[0,2], s=2, 
                axes_kwargs={'legend':False})
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='seq_run', ax=axs[1,0], s=2)
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='origin', ax=axs[1,1], s=2)
draw_embeddings(df, 'FA_diff1', 'FA_diff2', cat='condition', ax=axs[1,2], s=2)
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, '5000_no_unwanted_standard_FA_diff_subsampled.png'
    ), 
    dpi=500
)
# plt.show()

##

#########################################


##


######################################### Default, Harmony
a = AnnData(
    adata_original.layers['raw'], 
    obs=adata_original.obs.iloc[:,:11], 
    var=adata_original.var.iloc[:, [0,1]]
)
a.raw = a.copy()

# Seq run integration 5000 hvgs, auto pcs, Harmony, 30 k, 1 res, FA_diff.
_, reduced = standard_pp_recipe(a, n_HVGs=5000)
reduced = remove_unwanted(reduced)
d_pca = pca(reduced, layer='scaled', n_pcs=50)
reduced.obsm['scaled|original|X_pca'] = d_pca['X_pca']
compute_Harmony(reduced, categorical='seq_run')                                 # Harmony data int
compute_kNN(reduced, layer='scaled', int_method='Harmony', k=15)
reduced.obs['leiden_high_res'] = leiden_clustering(
    reduced.obsp['scaled|Harmony|X_corrected|NN_conn'], res=3.0
)
reduced.obs['leiden_high_res'] = pd.Categorical(reduced.obs['leiden_high_res']) # High res for embeddings
reduced.obs['leiden'] = leiden_clustering(
    reduced.obsp['scaled|Harmony|X_corrected|NN_conn'], res=.5               # Low res for clustering itself
)
reduced.obs['leiden'] = pd.Categorical(reduced.obs['leiden'])

# Embs
df = embeddings(reduced, rep='Harmony', paga_groups='leiden_high_res', umap_only=True)

# Cross
pd.crosstab(df['leiden'], df['seq_run'])
pd.crosstab(df['origin'], df['seq_run'])
pd.crosstab(df['sample'], df['seq_run'])

# fig, ax=plt.subplots()
# plot_heatmap(pd.crosstab(df['leiden'], df['seq_run']), ax=ax)
# plt.show()
# N.B. Good condition-seq-run corresponance

# Top clones
top_clones = df.groupby('GBC').size().loc[lambda x: x.apply(lambda x: x>50)].index
df['top_clones'] = df['GBC'].map(lambda x: x if x in top_clones else 'others')

# Plot
fig, axs = plt.subplots(2,3,figsize=(10,7))
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='leiden', ax=axs[0,0], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='top_clones', ax=axs[0,1], s=2,
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='sample', ax=axs[0,2], s=2, 
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='seq_run', ax=axs[1,0], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='origin', ax=axs[1,1], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='condition', ax=axs[1,2], s=2)
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, '5000_no_unwanted_standard_Harmony_UMAP_subsampled.png'
    ), 
    dpi=500
)
# plt.show()


##


######################################### no_cc, default
a = AnnData(
    adata_original.layers['raw'], 
    obs=adata_original.obs.iloc[:,:11], 
    var=adata_original.var.iloc[:, [0,1]]
)
a.raw = a.copy()

# Seq run integration 5000 hvgs, auto pcs, Harmony, 30 k, 1 res, FA_diff.
_, reduced = remove_cc_pp_recipe(a, n_HVGs=5000)
reduced = remove_unwanted(reduced)
make_folder(path_results, 'no_cc')
d_pca = pca(reduced, layer='scaled', n_pcs=50, auto=False, GSEA=True, path_viz=os.path.join(path_results, 'no_cc'))
reduced.obsm['scaled|original|X_pca'] = d_pca['X_pca']
compute_kNN(reduced, layer='scaled', k=15)
reduced.obs['leiden_high_res'] = leiden_clustering(
    reduced.obsp['scaled|original|X_pca|NN_conn'], res=3.0
)
reduced.obs['leiden_high_res'] = pd.Categorical(reduced.obs['leiden_high_res']) # High res for embeddings
reduced.obs['leiden'] = leiden_clustering(
    reduced.obsp['scaled|original|X_pca|NN_conn'], res=.6          # Low res for clustering itself
)
reduced.obs['leiden'] = pd.Categorical(reduced.obs['leiden'])

# Embs
df = embeddings(reduced, paga_groups='leiden_high_res', umap_only=True)

# Cross
pd.crosstab(df['leiden'], df['seq_run'])
pd.crosstab(df['origin'], df['seq_run'])
pd.crosstab(df['sample'], df['seq_run'])

# fig, ax=plt.subplots()
# plot_heatmap(pd.crosstab(df['leiden'], df['seq_run']), ax=ax)
# plt.show()
# N.B. Good condition-seq-run corresponance

# Top clones
top_clones = df.groupby('GBC').size().loc[lambda x: x.apply(lambda x: x>50)].index
df['top_clones'] = df['GBC'].map(lambda x: x if x in top_clones else 'others')

# Plot
fig, axs = plt.subplots(2,3,figsize=(10,7))
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='leiden', ax=axs[0,0], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='top_clones', ax=axs[0,1], s=2,
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='sample', ax=axs[0,2], s=2, 
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='seq_run', ax=axs[1,0], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='origin', ax=axs[1,1], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='condition', ax=axs[1,2], s=2)
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, '5000_no_unwanted_no_cc_UMAP_subsampled.png'
    ), 
    dpi=500
)
# plt.show()


##


######################################### no_cc, Harmony/Scanorama
a = AnnData(
    adata_original.layers['raw'], 
    obs=adata_original.obs.iloc[:,:11], 
    var=adata_original.var.iloc[:, [0,1]]
)
a.raw = a.copy()

# Seq run integration 5000 hvgs, auto pcs, Harmony, 30 k, 1 res, FA_diff.
_, reduced = remove_cc_pp_recipe(a, n_HVGs=5000)
reduced = remove_unwanted(reduced)
d_pca = pca(reduced, layer='scaled', n_pcs=50, auto=False)
reduced.obsm['scaled|original|X_pca'] = d_pca['X_pca']
compute_Scanorama(reduced, categorical='seq_run')                                 # Harmony data int
compute_kNN(reduced, layer='scaled', int_method='Scanorama', k=10)
reduced.obs['leiden_high_res'] = leiden_clustering(
    reduced.obsp['scaled|Scanorama|X_corrected|NN_conn'], res=2.0
)
reduced.obs['leiden_high_res'] = pd.Categorical(reduced.obs['leiden_high_res']) # High res for embeddings
reduced.obs['leiden'] = leiden_clustering(
    reduced.obsp['scaled|Scanorama|X_corrected|NN_conn'], res=.5             # Low res for clustering itself
)
reduced.obs['leiden'] = pd.Categorical(reduced.obs['leiden'])

# Embs
df = embeddings(reduced, rep='Scanorama', paga_groups='leiden_high_res', umap_only=True)

# Cross
pd.crosstab(df['leiden'], df['seq_run'])
pd.crosstab(df['origin'], df['seq_run'])
pd.crosstab(df['sample'], df['seq_run'])

# fig, ax=plt.subplots()
# plot_heatmap(pd.crosstab(df['leiden'], df['seq_run']), ax=ax)
# plt.show()
# N.B. Good condition-seq-run corresponance

# Top clones
top_clones = df.groupby('GBC').size().loc[lambda x: x.apply(lambda x: x>50)].index
df['top_clones'] = df['GBC'].map(lambda x: x if x in top_clones else 'others')

# Plot
fig, axs = plt.subplots(2,3,figsize=(10,7))
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='leiden', ax=axs[0,0], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='top_clones', ax=axs[0,1], s=2,
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='sample', ax=axs[0,2], s=2, 
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='seq_run', ax=axs[1,0], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='origin', ax=axs[1,1], s=2)
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='condition', ax=axs[1,2], s=2)
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, '5000_no_unwanted_no_cc_Scanorama_UMAP_subsampled.png'
    ), 
    dpi=500
)
# plt.show()


##


#########################################


##

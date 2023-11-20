"""
Cells states.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.gridspec import GridSpec
from plotting_utils._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula._utils import sanitize_neighbors
from BC_chemo_utils.utils import prep_df_for_dotplot
from BC_chemo_utils.tests import *
from BC_chemo_utils.plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'general_umaps')

# Read adata and format degs for dotplot
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
embs = (
    adata.obs
    .rename(columns={'cell_states':'cell_state'})
    .join(
    pd.DataFrame(
        adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names
    ))
)
df_markers = pd.read_csv(os.path.join(path_data, 'DEGs.csv'), index_col=0)

# Create colors
cell_state_colors = create_palette(embs, 'cell_state', 'tab20')
cell_state_colors['Undefined'] = '#E8E7E7'


##

# UMAP cell states + dotplot
fig = plt.figure(figsize=(15,5))
gs = GridSpec(1, 2, figure=fig, width_ratios=[2,2.5])

ax = fig.add_subplot(gs[0])
draw_embeddings(
    embs, cat='cell_state', ax=ax, title='Annotated cell states', 
    legend_kwargs={
        'ncols':1, 'colors':cell_state_colors, 
        'bbox_to_anchor':(1,1), 'loc':'upper left'
    },
)
# ax.axis('off')

ax = fig.add_subplot(gs[1])
df_ = prep_df_for_dotplot(df_markers)
df_['comparison'] = df_['comparison'].map(lambda x: x.split('_')[0])

sns.scatterplot(
    data=df_, y='comparison', x='gene', size='group_perc', hue='log2FC', 
    palette='mako', ax=ax, sizes=(1, 100)
)
format_ax(ax, title='Markers', xlabel='Top 3 marker genes', ylabel='Cell states', rotx=90)
ax.legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

fig.subplots_adjust(left=.05, right=.9, top=.9, bottom=.2, wspace=1.3)

# Save
fig.savefig(os.path.join(path_results, 'annotated_cell_states.png'), dpi=500)


##


# Enrichments

# 1. Origin
fig, axs = plt.subplots(1,2, figsize=(10, 4))

df_ = compute_enrichment(embs, 'cell_state', 'origin', 'PT')
sns.scatterplot(
    data=df_.sort_values('odds_ratio', ascending=False), 
    y='group', x='odds_ratio', size='perc_in_target', hue='FDR', 
    palette='Spectral', ax=axs[0], sizes=(1, 100)
)
format_ax(ax=axs[0], title='PT enrichment', xlabel='odds ratio', reduce_spines=True)
axs[0].legend().set_visible(False)

df_ = compute_enrichment(embs, 'cell_state', 'origin', 'lung')
sns.scatterplot(
    data=df_.sort_values('odds_ratio', ascending=False), 
    y='group', x='odds_ratio', size='perc_in_target', hue='FDR', 
    palette='Spectral', ax=axs[1], sizes=(1, 100)
)
format_ax(ax=axs[1], title='Lung enrichment', xlabel='odds ratio', reduce_spines=True)
axs[1].legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'origin_enrichment.png'), dpi=300)


##


# 1. Condition
fig = plt.figure(figsize=(15, 7))

for i,t in enumerate(embs['condition'].cat.categories):

    df_ = compute_enrichment(embs, 'cell_state', 'condition', t)
    ax = fig.add_subplot(2,3,i+1)
    sns.scatterplot(
        data=df_.sort_values('odds_ratio', ascending=False), 
        y='group', x='odds_ratio', size='perc_in_target', hue='FDR', 
        palette='Spectral', ax=ax, sizes=(1, 100)
    )
    max_ = df_['odds_ratio'].max()
    delta = max_ + max_ * 0.01
    ax.set_xlim((-delta/5, max_ + delta))
    format_ax(ax=ax, title=f'{t} enrichment', xlabel='odds ratio', reduce_spines=True)
    if i+1 != 6:
        ax.legend().set_visible(False)
    else:
        ax.legend(loc='center left', bbox_to_anchor=(2,1.2), frameon=False)

fig.subplots_adjust(left=.2, right=.7, top=.9, bottom=.1, wspace=2.3, hspace=.4)
fig.savefig(os.path.join(path_results, 'condition_enrichment.png'), dpi=300)


##


# PAGA
adata = sanitize_neighbors(adata)
sc.tl.paga(adata, groups='cell_states', neighbors_key='nn')

# Network
fig, ax = plt.subplots(figsize=(7,7))
sc.pl.paga(adata, frameon=False, ax=ax)
sc.pl.paga(adata, frameon=False, ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'paga.png'), dpi=1000)

# Heatmap
cats = adata.obs['cell_states'].cat.categories
df_ = pd.DataFrame(adata.uns['paga']['connectivities'].A, index=cats, columns=cats)
g = plot_clustermap(df_, figsize=(7,5), cb_label='Connectivities', title='PAGA cell states')
g.fig.savefig(os.path.join(path_results, 'paga_heat.png'), dpi=300)


##A


# CC

fig = plt.figure(figsize=(14,6.8))

gs = GridSpec(2, 6, figure=fig, height_ratios=[1,1.5])

ax = fig.add_subplot(gs[0,1:-1])
pairs = [
    ['PT, untreated', 'PT, treated'],
    ['lung, untreated', 'lung, PT-treated'],
    ['lung, untreated', 'lung, lung-treated'],
    ['lung, lung-treated', 'lung, double-treated'],
    ['lung, PT-treated', 'lung, double-treated'],
]
violin(embs, 'condition', 'cycling', ax=ax, c='darkgrey', with_stats=True, pairs=pairs)
format_ax(ax, title='Cell cycle signatures scores', 
          xticks=embs['condition'].cat.categories, ylabel='Score')
ax.spines[['left', 'right', 'top']].set_visible(False)

ax = fig.add_subplot(gs[1,:2])
draw_embeddings(embs, cont='s_seurat', ax=ax, title='cycling', cbar_kwargs={'palette':'mako'}, s=1)
ax.axis('off')
ax = fig.add_subplot(gs[1,2:4])
draw_embeddings(embs, cont='G1/S', ax=ax, title='G1/S', cbar_kwargs={'palette':'mako'}, s=1)
ax.axis('off')
ax = fig.add_subplot(gs[1,4:])
draw_embeddings(embs, cont='G2/M', ax=ax, title='G2/M', cbar_kwargs={'palette':'mako'}, s=1)
ax.axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cc.png'), dpi=400)


##
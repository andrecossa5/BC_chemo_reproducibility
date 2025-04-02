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
path_main = '/Users/ieo7295/Desktop/BC_chemo_reproducibility/'
path_adata = '/Users/ieo7295/Desktop/tests/cell/data/default'
path_data = os.path.join(path_main, 'data', 'CTCs')
path_results = os.path.join(path_main, 'results', 'CTCs', 'general_umaps')

# Read adata and format degs for dotplot
adata = sc.read(os.path.join(path_adata, 'clustered.h5ad'))
embs = (
    adata.obs
    .rename(columns={'final_cell_states':'cell_state'})
    .join(
    pd.DataFrame(
        adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names
    ))
)

# Create colors
cell_state_colors = create_palette(embs, 'final_cell_state', 'tab20')
cell_state_colors['Undefined'] = '#E8E7E7'


# Enrichments

# 1. Origin
fig, axs = plt.subplots(1,3, figsize=(16, 4))

df_ = compute_enrichment(embs, 'final_cell_state', 'origin', 'PT')
sns.scatterplot(
    data=df_.sort_values('odds_ratio', ascending=False), 
    y='group', x='odds_ratio', size='perc_in_target', hue='FDR', 
    palette='Spectral', ax=axs[0], sizes=(1, 100)
)
format_ax(ax=axs[0], title='PT enrichment', xlabel='odds ratio', reduce_spines=True)
axs[0].legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

df_ = compute_enrichment(embs, 'final_cell_state', 'origin', 'lung')
sns.scatterplot(
    data=df_.sort_values('odds_ratio', ascending=False), 
    y='group', x='odds_ratio', size='perc_in_target', hue='FDR', 
    palette='Spectral', ax=axs[1], sizes=(1, 100)
)
format_ax(ax=axs[1], title='Lung enrichment', xlabel='odds ratio', reduce_spines=True)
axs[1].legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

df_ = compute_enrichment(embs, 'final_cell_state', 'origin', 'CTC')
sns.scatterplot(
    data=df_.sort_values('odds_ratio', ascending=False), 
    y='group', x='odds_ratio', size='perc_in_target', hue='FDR', 
    palette='Spectral', ax=axs[2], sizes=(1, 100)
)
format_ax(ax=axs[2], title='CTC enrichment', xlabel='odds ratio', reduce_spines=True)
axs[2].legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'origin_enrichment.png'), dpi=300)


##


# PAGA
adata = sanitize_neighbors(adata)
sc.tl.paga(adata, groups='final_cell_state', neighbors_key='nn')

# Network
fig, ax = plt.subplots(figsize=(7,7))
sc.pl.paga(adata, frameon=False, ax=ax)
sc.pl.paga(adata, frameon=False, ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'paga.png'), dpi=1000)

# Heatmap
cats = adata.obs['final_cell_state'].cat.categories
df_ = pd.DataFrame(adata.uns['paga']['connectivities'].A, index=cats, columns=cats)
g = plot_clustermap(df_, figsize=(7,5), cb_label='Connectivities', title='PAGA cell states')
g.fig.savefig(os.path.join(path_results, 'paga_heat.png'), dpi=300)


##


# CC

fig = plt.figure(figsize=(14,6.8))

gs = GridSpec(2, 6, figure=fig, height_ratios=[1,1.5])

ax = fig.add_subplot(gs[0,1:-1])
pairs = [
    ['CTC', 'PT'],
    ['CTC', 'lung'],
    ['PT', 'lung'],
]

def violin(df, x, y, by=None, c=None, a=1, l=None, ax=None, with_stats=False, order=None, pairs=None):
    """
    Base violinplot.
    """
    params = {   
        'boxprops' : {'edgecolor': 'white', 'linewidth': 0.5}, 
        'medianprops': {"color": "white", "linewidth": 1.2},
        'whiskerprops':{"color": "black", "linewidth": 1}
    }
    
    if isinstance(c, str):
        ax = sns.violinplot(data=df, x=x, y=y, color=c, ax=ax, saturation=0.7, order=order, **params) 
        ax.set(xlabel='', ylabel='')
        ax.set_xticklabels(np.arange(df[x].unique().size))

    elif isinstance(c, dict) and by is None:
        ax = sns.violinplot(data=df, x=x, y=y, palette=c.values(), ax=ax, saturation=0.7, order=order, **params)
        ax.set(xlabel='', ylabel='') 
        ax.set_xticklabels(np.arange(df[x].unique().size))
            
    elif isinstance(c, dict) and by is not None:
        ax = sns.violinplot(data=df, x=x, y=y, palette=c.values(), hue=by, 
            ax=ax, saturation=0.7, **params)
        ax.legend([], [], frameon=False)
        ax.set(xlabel='', ylabel='')
        ax.set_xticklabels(np.arange(df[x].unique().size))

    else:
        raise ValueError(f'{by} categories do not match provided colors keys')

    if with_stats:
        add_wilcox(df, x, y, pairs, ax, order=None)

    return ax

violin(adata.obs, 'origin', 'cycling', ax=ax, c='darkgrey', with_stats=True, pairs=pairs)
format_ax(ax, title='Cell cycle signatures scores', 
          xticks=adata.obs['origin'].cat.categories, ylabel='Score')
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
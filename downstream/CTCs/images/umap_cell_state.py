"""
Basic UMAPs on transcriptional space.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
from Cellula.plotting._plotting import *
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *
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
path_main = '/Users/ieo7295/Desktop/tests/cell'
path_adata = os.path.join(path_main, 'data', 'default')
path_results = '/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/CTCs/general_umaps'

##

# Read adata and format degs for dotplot
adata = sc.read(os.path.join(path_adata, 'clustered.h5ad'))
embs = (
    adata.obs
    .rename(columns={'final_cell_state':'cell_state'})
    .join(
    pd.DataFrame(
        adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names
    ))
)
df_markers = pd.read_csv(os.path.join(path_results, 'DEGs.csv'), index_col=0)
# Create colors
cell_state_colors = create_palette(embs, 'cell_state', 'tab20')
cell_state_colors['Undefined'] = '#E8E7E7'

## UMAP cell_state + dotplot
fig = plt.figure(figsize=(20,8))
gs = GridSpec(1, 2, figure=fig, width_ratios=[2,2.5])

ax = fig.add_subplot(gs[0])
draw_embeddings(
    embs, cat='cell_state', ax=ax, title='Annotated cell states', 
    legend_kwargs={
        'ncols':1, 'colors':cell_state_colors, 
        'bbox_to_anchor':(1,1), 'loc':'upper left'
    },
)

ax = fig.add_subplot(gs[1])
df_ = prep_df_for_dotplot(df_markers)
df_['comparison'] = df_['comparison'].map(lambda x: x.split('_')[0][:-1])

sns.scatterplot(
    data=df_, y='comparison', x='gene', size='group_perc', hue='log2FC', 
    palette='mako', ax=ax, sizes=(1, 100)
)
format_ax(ax, title='Markers', xlabel='Top 3 marker genes', ylabel='Cell states', rotx=90)
ax.legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

fig.subplots_adjust(left=.05, right=.9, top=.9, bottom=.2, wspace=1.3)
plt.show()
# Save
fig.savefig(os.path.join(path_results, 'annotated_cell_states.png'), dpi=500)


#PAGA
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



# CC

fig = plt.figure(figsize=(14,6.8))

gs = GridSpec(2, 6, figure=fig, height_ratios=[1,1.5])

ax = fig.add_subplot(gs[0,1:-1])
pairs = [
    ['PT', 'lung'],
    ['CTC', 'lung'],
    ['CTC', 'PT'],
]


def violin(df, x, y, by=None, c=None, a=1, l=None, ax=None, with_stats=False, order=None, pairs=None):
    """
    Base violinplot, updated for seaborn >= 0.12 and matplotlib >= 3.8.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    plot_kwargs = {
        'data': df,
        'x': x,
        'y': y,
        'ax': ax,
        'saturation': a,
        'order': order,
        'linewidth': 1
    }

    # Handle color/palette logic
    if isinstance(c, str):
        plot_kwargs['color'] = c
    elif isinstance(c, dict):
        plot_kwargs['palette'] = [c[k] for k in df[x].unique() if k in c]

    # Plot
    if by is None:
        sns.violinplot(**plot_kwargs)
    else:
        plot_kwargs['hue'] = by
        sns.violinplot(**plot_kwargs)
        ax.legend([], [], frameon=False)

    ax.set(xlabel='', ylabel='')

    # Optional stats annotation
    if with_stats and pairs:
        try:
            add_wilcox(df, x, y, pairs, ax, order=order)
        except Exception as e:
            print(f"[Warning] Failed to add statistical annotations: {e}")

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

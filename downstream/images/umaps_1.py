"""
Clonal summary stats plots.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from plotting_utils._utils import *
from plotting_utils._colors import *
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
from Cellula.plotting._plotting import draw_embeddings, faceted_draw_embedding
from matplotlib.gridspec import GridSpec
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'general_umaps')

# Read cells meta and format
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
embs = adata.obs.join(
    pd.DataFrame(
        adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs_names
    )
)
top_clones = embs.groupby('GBC').size().loc[lambda x:x>100].index
embs['GBC_top'] = embs['GBC'].map(lambda x: x if x in top_clones else 'others')


##


# Read clone colors and modify it
clone_colors = pd.read_csv(os.path.join(path_data, 'clones_colors_sc.csv'), index_col=0)
clone_colors = clone_colors['color'].to_dict()
clone_colors = { k:clone_colors[k] for k in embs['GBC_top'] if k != 'others'}
clone_colors['others'] = '#E8E7E7'

# Condition colors
condition_colors = create_palette(embs, 'condition', col_list=ten_godisnot)


##


# 1. UMAP1: clones, category, origin, nUMIs, mito_perc, cycling
fig = plt.figure(figsize=(15,5))
gs = GridSpec(2, 4, figure=fig, width_ratios=[2,1,1,1], height_ratios=[1,1])

ax = fig.add_subplot(gs[:,0])
draw_embeddings(
    embs, cat='GBC_top', ax=ax, title=f'Top clones (n={len(clone_colors)-1})', 
    legend_kwargs={'color':clone_colors},
    axes_kwargs={'legend':False}
)
# ax.axis('off')

ax = fig.add_subplot(gs[0,1])
draw_embeddings(
    embs, cat='origin', ax=ax, title='Origin',
    legend_kwargs={'bbox_to_anchor':(.9,1), 'loc':'upper left'},    
)
ax.axis('off')

ax = fig.add_subplot(gs[0,2])
draw_embeddings(
    embs, cat='condition', ax=ax, title='Condition', 
    legend_kwargs={
        'ncols':2, 'colors':condition_colors, 
        'bbox_to_anchor':(.9,1), 'loc':'upper left'
    },
)
ax.axis('off')

ax = fig.add_subplot(gs[1,1])
draw_embeddings(embs, cont='nUMIs', ax=ax, title='nUMIs', cbar_kwargs={'palette':'mako'})
ax.axis('off')

ax = fig.add_subplot(gs[1,2])
draw_embeddings(
    embs, cont='mito_perc', ax=ax, title='% mito',
    cbar_kwargs={'vmin':.05, 'vmax':.1, 'palette':'mako'}
)
ax.axis('off')

ax = fig.add_subplot(gs[1,3])
draw_embeddings(embs, cont='cycling', ax=ax, title='Cycling signature', cbar_kwargs={'palette':'mako'})
ax.axis('off')

fig.subplots_adjust(top=.9, bottom=.1, left=.1, right=.9, wspace=.3, hspace=.2)
fig.savefig(os.path.join(path_results, 'UMAP_1.png'), dpi=500)


##



colors = create_palette(embs, 'sample', sc.pl.palettes.default_102)

fig, ax = plt.subplots(figsize=(10.5,7))
draw_embeddings(
    embs, cat='seq_run', ax=ax, title='', s=3,
    legend_kwargs={
        'ncols':1,
        'bbox_to_anchor':(.9,1), 'loc':'upper left'
    },
)
ax.axis('off')
fig.tight_layout()
plt.show()
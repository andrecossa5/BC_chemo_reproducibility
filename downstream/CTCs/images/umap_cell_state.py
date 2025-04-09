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

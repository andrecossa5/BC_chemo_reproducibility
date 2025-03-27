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
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/ieo7295/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'CTCs')
path_results = os.path.join(path_main, 'results', 'CTCs', 'general_umaps')


# Read cells meta and format
adata =sc.read_h5ad(os.path.join(path_data, 'clustered.h5ad'))

##


# Extract embeddings
df_umap= pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1','UMAP2'], index=adata.obs.index)
df_u = pd.concat([adata.obs, df_umap], axis=1)

#Cell colors
cell_state_colors = create_palette(df_u, 'final_cell_state', col_list=sc.pl.palettes.default_20)
cell_state_colors['Undefined'] = 'grey'
df_u = df_u.rename(columns={'final_cell_state':'Cell states'})


# Viz
fig, ax = plt.subplots(figsize=(6.5,4.5))
draw_embeddings(
    df_u, cat='Cell states', ax=ax, title=f'Cell states', 
    legend_kwargs={
        'colors':cell_state_colors,
        'bbox_to_anchor' : (1,.5),
        'loc' : 'center left',
    },
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'UMAP_cell_states.png'), dpi=350)

##



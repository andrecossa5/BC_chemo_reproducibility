"""
Basic UMAPs on transcriptional space.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'CTCs')
path_results = os.path.join(path_main, 'results', 'CTCs', 'general_umaps')


# Read cells meta and format
embs = pd.read_csv(os.path.join(path_data, 'full_embs.csv'), index_col=0)


##


# Cell state colors
cell_state_colors = create_palette(embs, 'final_cell_state', col_list=sc.pl.palettes.default_20)
cell_state_colors['Undefined'] = 'grey'
embs = embs.rename(columns={'final_cell_state':'Cell state'})


##


# Viz
fig, ax = plt.subplots(figsize=(6.5,4.5))
draw_embeddings(
    embs, cat='Cell state', ax=ax, title=f'Cell states', 
    legend_kwargs={
        'colors':cell_state_colors,
        'bbox_to_anchor' : (1,.5),
        'loc' : 'center left',
    },
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'UMAP_cell_states.png'), dpi=500)


##
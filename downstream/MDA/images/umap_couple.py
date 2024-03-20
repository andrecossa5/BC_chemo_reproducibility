"""
Couple umap
"""

import os
import pickle
import pandas as pd
from plotting_utils._utils import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'general_umaps')

# Read cells meta and umap
embs = pd.read_csv(os.path.join(path_data, 'full_embs.csv'), index_col=0)


##


# Read clone colors and modify it
with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'rb') as f:
    clone_colors = pickle.load(f)

# Create the color map using Seaborn
sns.color_palette("Spectral_r")
origin_colors = { 
    "PT" : sns.color_palette("Spectral_r")[0], 
    "lung" : sns.color_palette("Spectral_r")[-1] 
}

##


# UMAP all
fig, ax = plt.subplots(figsize=(4.5,5))
draw_embeddings(
    embs, 
    cat='origin', 
    ax=ax,
    legend_kwargs={
        'bbox_to_anchor':(1,1), 
        'loc':'upper right', 
        'ncols':1,
        'label':'',
        'artists_size' : 11,
        'label_size' : 15,
        'ticks_size' : 13,
        'colors' : origin_colors
    }
)
ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'UMAP_all_origin.pdf'), dpi=500)


##


# UMAP NT_NT_2
fig, ax = plt.subplots(figsize=(4.5,5))
draw_embeddings(
    embs.query('dataset=="NT_NT_2"'), 
    cat='origin', 
    ax=ax,
    axes_kwargs={'legend':False},
    legend_kwargs={
        'bbox_to_anchor':(1,1), 
        'loc':'upper right', 
        'ncols':1,
        'label':'',
        'label_size' : 15,
        'ticks_size' : 13,
        'colors' : origin_colors
    }
)
ax.axis('off')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'UMAP_NT_NT2_origin.pdf'), dpi=500)


##

##
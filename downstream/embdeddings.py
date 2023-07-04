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

# Read 
clustered = sc.read(os.path.join(path_data, 'subsampled.h5ad'))

# Manual clustering
# for r in [.2,.3,.4,.5,.6,.7,.8,.9,1]:
#     leiden(clustered, obsp_key='NN_conn', obs_key=f'leiden_{r}', res=r)
#     df = clustered.obs.join(
#         pd.DataFrame(
#             clustered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=clustered.obs_names
#         )
#     )
#     fig, ax = plt.subplots(figsize=(5,5))
#     draw_embeddings(df, 'UMAP1', 'UMAP2', cat=f'leiden_{r}', ax=ax)
#     ax.axis('off')
#     fig.tight_layout()
#     fig.savefig(
#         os.path.join(
#             path_results, 'manual_clustering', f'leiden_{r}.png'
#         ), 
#         dpi=500
#     )

# Create df for plotting
df = clustered.obs.join(
    pd.DataFrame(clustered.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=clustered.obs_names)
)

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

# Cat 
fig, axs = plt.subplots(2,3,figsize=(10,7))
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='leiden', ax=axs[0,0])
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='top_clones', ax=axs[0,1],
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='sample', ax=axs[0,2],
                axes_kwargs={'legend':False})
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='seq_run', ax=axs[1,0])
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='origin', ax=axs[1,1])
draw_embeddings(df, 'UMAP1', 'UMAP2', cat='condition', ax=axs[1,2])
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, 'Definitive_embs_cat.png'
    ), 
    dpi=500
)

# Cont
fig, axs = plt.subplots(1,4,figsize=(14.5,3.5))
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='cycling', ax=axs[0],
                cbar_kwargs={'vmin':0, 'vmax':.5, 'layout':( (1.05,.25,.03,.5), 'right' )})
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='G1/S', ax=axs[1],
                cbar_kwargs={'vmin':0, 'vmax':.2, 'layout':( (1.05,.25,.03,.5), 'right' )})
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='G2/M', ax=axs[2],
               cbar_kwargs={'vmin':0, 'vmax':.5, 'layout':( (1.05,.25,.03,.5), 'right' )})
draw_embeddings(df, 'UMAP1', 'UMAP2', cont='ribo_genes', ax=axs[3],
                cbar_kwargs={'vmin':.7, 'vmax':4, 'layout':( (1.05,.25,.03,.5), 'right' )})
fig.tight_layout()
fig.savefig(
    os.path.join(
        path_results, 'Definitive_embs_cont.png'
    ), 
    dpi=500
)


##


# Custom UMAP, AC_AC 4 dominant 
fig, axs = plt.subplots(4,2,figsize=(4.5, 9))

# AC_AC_1
df_ = (
    df
    .query('dataset == "AC_AC_1"')
    .assign(clone=lambda x: np.where(x['GBC']=='CGAGGGGATGGACTTCCG', 'clone', 'rest'))
)
draw_embeddings(df_.query('origin == "PT"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[0,0], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#CC0000', 'rest':'#C0C0C0'}},
    title='Clone CGAGGG (PT)'      
)
draw_embeddings(df_.query('origin == "lung"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[0,1], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#CC0000', 'rest':'#C0C0C0'}},
    title='Clone CGAGGG (lung)'
)

# AC_AC_2
df_ = (
    df
    .query('dataset == "AC_AC_2"')
    .assign(clone=lambda x: np.where(x['GBC']=='AACTCGACGCCTTATCAG', 'clone', 'rest'))
)
draw_embeddings(df_.query('origin == "PT"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[1,0], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#00CC66', 'rest':'#C0C0C0'}},
    title='Clone AACTCG (PT)'      
)
draw_embeddings(df_.query('origin == "lung"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[1,1], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#00CC66', 'rest':'#C0C0C0'}},
    title='Clone AACTCG (lung)'
)

# AC_AC_2
df_ = (
    df
    .query('dataset == "AC_AC_2"')
    .assign(clone=lambda x: np.where(x['GBC']=='GAAAGCGTAACGCGTCAG', 'clone', 'rest'))
)
draw_embeddings(df_.query('origin == "PT"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[2,0], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#CC6600', 'rest':'#C0C0C0'}},
    title='Clone GAAAGC (PT)'      
)
draw_embeddings(df_.query('origin == "lung"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[2,1], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#CC6600', 'rest':'#C0C0C0'}},
    title='Clone GAAAGC (lung)'
)

# AC_AC_3
df_ = (
    df
    .query('dataset == "AC_AC_3"')
    .assign(clone=lambda x: np.where(x['GBC']=='GTTCTTAGGTGTCCAGGT', 'clone', 'rest'))
)
draw_embeddings(df_.query('origin == "PT"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[3,0], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#009999', 'rest':'#C0C0C0'}},
    title='Clone GTTCTT (PT)'      
)
draw_embeddings(df_.query('origin == "lung"'),
    'UMAP1', 'UMAP2', cat='clone', ax=axs[3,1], 
    axes_kwargs={'legend':False},
    legend_kwargs={'colors': {'clone':'#009999', 'rest':'#C0C0C0'}},
    title='Clone GTTCTT (lung)'
)

# Save
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'UMAP_cloni_ACAC.png'), dpi=500)


##



"""
Default, MDA231 chemo, 1 Feb 2023.
"""

# Code
import sys
import os
from Cellula._utils import *
from Cellula.preprocessing._pp import * 
from Cellula.preprocessing._neighbors import * 
from Cellula.preprocessing._embeddings import *
from Cellula.preprocessing._integration import * 
from Cellula.clustering._clustering import * 
from Cellula.plotting._plotting import *
from Cellula.plotting._colors import *
matplotlib.use('MacOSX')


# Paths
path_main = '/Users/IEO5505/Desktop/MDA_chemo_repro/single_cell/'
make_folder(path_main + 'results_and_plots/viz/', 'categoricals')
make_folder(path_main + 'results_and_plots/report/', 'categoricals')

path_report = path_main + 'results_and_plots/report/categoricals/'
path_viz = path_main + 'results_and_plots/viz/categoricals/'


##


# Load data 

# Data
adata = sc.read(path_main + 'data/lognorm.h5ad') #clustered.h5ad
#embs = pd.read_csv(path_main + 'data/embeddings.csv', index_col=0)

# NB: Remove <-10-cells-per-sample-clones' cells
adata = adata[adata.obs['GBC_10_per_sample'] != "others", :].copy()

## Pp
t = Timer()

t.start()
adata = red(adata) # add n_hvgs=2000 kwarg + docstring
adata = scale(adata)
adata = pca(adata)
adata = compute_kNN(adata)
df = embeddings(adata, layer='scaled')
t.stop()

## Clustering
t.start()
for r in np.round(np.linspace(0.5,1.5,10),2):
    adata.obs[f'Leiden_{r}'] = leiden_clustering(
        adata.obsp['scaled|original|X_pca|15_NN_30_comp_conn'],
    )
    print(adata.obs[f'Leiden_{r}'])
t.stop()

# Join 
df = adata.obs.join(df)

# Top clones
top = df.groupby('GBC_10_per_sample').size(
    ).sort_values(ascending=False)[:50].index
df['top_50'] = [ x if x in top else 'others' for x in df['GBC_10_per_sample'] ]


##


################################################################

# Categorical plots

# Fig, umap categorical
# for x in df.columns:
#     if df[x].dtype in ['category', 'object'] and df[x].dtype != 'GBC':
#         print(x)
#         fig, ax = plt.subplots(figsize=(8.5,8.5))
#         draw_embeddings(df, cat=x, ax=ax, s=1)
#         #ax.xaxis.set_tick_params(labelsize=15)
#         #ax.yaxis.set_tick_params(labelsize=15)
#         fig.savefig(path_viz + f'{x}_umap.pdf')


##

fig, ax = plt.subplots(figsize=(8.5, 6))

draw_embeddings(df, cat='top_50', ax=ax, s=0.1, title='Prova clones',
    legend_kwargs={
        'bbox_to_anchor' : (1.02,.5),
        'loc' : 'center left', 
        'label_size' : 7,
        'ticks_size' : 5,
        'label' : 'Clones',
        'colors' : create_palette(df, 'top_50', sc.pl.palettes.godsnot_102),
        'only_top' : 'all',
        'ncols' : 2
    },
    axes_kwargs={
        'title_size' : 13,
        'xlabel_size' : 12,
        'ylabel_size' : 12,
        'legend' : True
    },
    query='origin == "PT"'
)

# draw_embeddings(df, cont='cycling', ax=axs[1], s=0.1, title='Cycl',
#     cbar_kwargs={
#         'pos' : 'outside',
#         'orientation' : 'v',
#         'label_size' : 10,
#         'ticks_size' : 8,
#         'color' : 'viridis',
#     },
#     axes_kwargs={
#         'title_size' : 20,
#         'xlabel_size' : 2,
#         'ylabel_size' : 25,
#     },
#     query='origin == "PT"'
# )

plt.subplots_adjust(right=0.6, top=0.8, bottom=0.15, left=0.15)
# fig.tight_layout()

fig.savefig(path_viz + 'top_50_clones_PT.pdf')


# s_top = df.loc[df['GBC'].isin(top)].groupby(['sample', 'GBC']).size()
# s_top = s[s!=0].to_frame('size').reset_index()
# 
# 
# 
# freq_top = pd.merge(
#     s_top.reset_index(),
#     df.loc[:, ['sample', 'origin', 'dataset']].reset_index(drop=True).drop_duplicates(),
#     how='left',
#     on='sample'
# ).loc[:, ['GBC', 'origin', 'dataset', 'size']]
# 
# freq_top.to_excel(path_report + 'top10.xlsx')


##





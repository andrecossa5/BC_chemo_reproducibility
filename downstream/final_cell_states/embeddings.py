"""
Visualize embeddings.
"""

import os
from Cellula.plotting._colors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'final_cell_states')


##


# Data
adata = sc.read(os.path.join(path_data, 'subsample.h5ad'))
df = pd.read_csv(os.path.join(path_data, 'full_embs.csv'), index_col=0)
cons_clust = pd.read_csv(os.path.join(path_data, '50_NN_0.94_cons_clusters.csv'), index_col=0)
df['consensus'] = cons_clust['0']

# Summary QC
(
    df.loc[:, ['nUMIs', 'detected_genes', 'mito_perc', 'cell_state_lowres']]
    .groupby('cell_state_lowres')
    .agg('median')
)


##


# Viz embeddings
make_folder(path_results, 'visualization', overwrite=True)

cats = ['consensus', 'leiden', 'seq_run', 'origin']
conts = ['nUMIs', 'mito_perc', 'detected_genes', 
        's_seurat', 'g2m_seurat', 'cycle_diff', 'cycling']

# Single covariates
for x in cats:
    fig, ax = plt.subplots(figsize=(5,5))
    draw_embeddings(df, cat=x, ax=ax, axes_kwargs={'legend':False})
    ax.axis('off')
    add_labels_on_loc(df, 'UMAP1', 'UMAP2', x, ax, s=7)
    fig.tight_layout()
    fig.savefig(
        os.path.join(path_results, 'visualization', f'{x}.png'),
        dpi=300
    )
for x in conts:
    fig, ax = plt.subplots(figsize=(5.5,5))
    draw_embeddings(df, cont=x, ax=ax)
    ax.axis('off')
    fig.tight_layout()
    fig.savefig(
        os.path.join(path_results, 'visualization', f'{x}.png'),
        dpi=300
    )
  
# Faceted by condition
for x in cats:
    fig = faceted_draw_embedding(
        df, cat=x, facet='condition', figsize=(13,4), n_cols=4,
        lables_on_loc=True, axis=False, legend=False
    )
    fig.savefig(
        os.path.join(path_results, 'visualization', f'{x}_by_condition.png'),
        dpi=300
    )


##


# Manual annot
fig, ax = plt.subplots(figsize=(9,5))
colors = create_palette(df, 'cell_state_highres', sc.pl.palettes.vega_20)
draw_embeddings(
    df, cat='cell_state_highres', ax=ax, title='Manual annotation, high resolution',
    legend_kwargs={'colors':colors, 'loc':'upper left', 'bbox_to_anchor':(1,1)},
)
fig.subplots_adjust(left=.1, bottom=.1, top=.9, right=.6)
ax.axis('off')
fig.tight_layout()
fig.savefig(
    os.path.join(path_results, 'visualization', 'manual_high_resolution.png'), dpi=300
)

fig, ax = plt.subplots(figsize=(9,5))
bb_plot(
    df, 'condition', 'cell_state_highres', colors=colors, 
    ax=ax, bbox_to_anchor=(1,1), loc='upper left', 
    label_size=10, ticks_size=8, artists_size=8
)
fig.subplots_adjust(left=.1, bottom=.1, top=.9, right=.6)
fig.savefig(
    os.path.join(path_results, 'visualization', 'condition_bbplot.png'), dpi=300
)

fig, ax = plt.subplots(figsize=(9,5))
bb_plot(
    df, 'origin', 'cell_state_highres', colors=colors, 
    ax=ax, bbox_to_anchor=(1,1), loc='upper left', 
    label_size=10, ticks_size=8, artists_size=8
)
fig.subplots_adjust(left=.1, bottom=.1, top=.9, right=.6)
fig.savefig(
    os.path.join(path_results, 'visualization', 'origin_bbplot.png'), dpi=300
)

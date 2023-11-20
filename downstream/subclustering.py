#!/usr/bin/python

import os
import pandas as pd
import scanpy as sc
from anndata import AnnData
from Cellula.preprocessing._integration import harmonize, sc
from Cellula.plotting._plotting import *


# ...

a = adata[~adata.obs['dataset'].str.startswith('AC_NT'),:].copy()


# for dataset in adata.obs['dataset'].unique():

a = AnnData(X=a.layers['raw'], obs=a.obs, var=a.var)


a.X[:10,:10]

sc.pp.normalize_total(a, target_sum=10**4)
sc.pp.highly_variable_genes(a, n_top_genes=3000, flavor='cell_ranger')
sc.pp.scale(a)
sc.pp.regress_out(a, keys=['cycle_diff'], n_jobs=8, copy=True)
sc.pp.scale(a)
sc.pp.pca(a, n_comps=10)
# a.obsm['X_harmony'] = harmonize(a.obsm['X_pca'], a.obs, 'sample')
sc.pp.neighbors(a, use_rep='X_pca', n_neighbors=15)
sc.tl.leiden(a, resolution=.3)
sc.tl.paga(a, groups='leiden')
sc.pl.paga(a, plot=False)
sc.tl.umap(a, init_pos='paga')


##

_embs = (
    a.obs
    .join(
    pd.DataFrame(
        a.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=a.obs_names
    ))
)
for c in _embs.columns:
    if pd.api.types.is_categorical_dtype(_embs[c]):
        _embs[c] = _embs[c].astype('str')

    ##


fig, axs = plt.subplots(1,2,figsize=(9,4))
draw_embeddings(_embs, cat='leiden', ax=axs[0], title='Leiden', s=3,
                legend_kwargs={'loc':'upper left'})
axs[0].axis('off')
draw_embeddings(_embs, cat='sample', ax=axs[1], title='Sample', s=3,
                legend_kwargs={'loc':'upper left'})
axs[1].axis('off')
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'NO_AC_NT.png'), dpi=400)


# draw_embeddings(_embs, cat='cell_states', ax=axs[1], title='Cell_state', s=1, 
#                 legend_kwargs={'loc':'upper left'})
# axs[1].axis('off')
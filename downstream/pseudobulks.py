"""
Pseudobulk analysis.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import pegasus as pg
import matplotlib
import matplotlib.pyplot as plt

from itertools import product
from fbpca import pca
from pegasusio import MultimodalData
from Cellula.plotting._colors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


## Utils
def simple_pca(df):
    """
    Simple PCA for speudobulk data
    """
    
    np.random.seed(1234)
    k = 10 if 10 < df.shape[0] else df.shape[0]

    from sklearn.preprocessing import scale
    X_rescaled = scale(df.values)

    X_pca, sv, loads = pca(
        X_rescaled,
        k=k, 
        raw=True
    )
    loads = loads.T
    sqvars = sv**2
    var_ratios = sqvars / np.sum(sqvars)
    cum_sum_eigenvalues = np.cumsum(var_ratios)

    loads = pd.DataFrame(
        loads[:,:5],
        index=df.columns,
        columns=[ f'PC{i+1}' for i in range(5) ]
    )
    df_pcs = pd.DataFrame(X_pca, index=df.index, columns=[ f'PC{i+1}' for i in range(k) ])

    return df_pcs, loads, cum_sum_eigenvalues


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'pseudobulk')


#

# Read data 
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))

# Remove AC_NT branch
adata = adata[adata.obs.query('condition != "AC_NT"').index,:]

# Sample
data = pg.pseudobulk(MultimodalData(adata), 'sample', mat_key='raw')
counts = pd.DataFrame(data.matrices['counts'], index=data.obs.index, columns=data.var.index)
lognorm = counts.apply(lambda x: x*1000000 / x.sum(), axis=1) # CPM
X_pca, loads, cum_sum_eigenvalues = simple_pca(lognorm)

# Viz
df_ = (
    X_pca
    .join(
        adata.obs.loc[:,['condition', 'origin', 'sample', 'dataset']]
        .drop_duplicates().set_index('sample')
    )
    .reset_index()
    .rename(columns={'barcodekey':'sample'})
)
colors = {
    'condition' : create_palette(df_, 'condition', ten_godisnot),
    'sample' : create_palette(df_, 'sample', sc.pl.palettes.default_28),
    'origin' : create_palette(df_, 'origin', ten_godisnot[::-4]),
    'dataset' : create_palette(df_, 'dataset', ten_godisnot[::-1])
}

# Viz
pcs = df_.columns[df_.columns.str.contains('PC')][:5]
for x, y in product(pcs, pcs):

    if x != y:
        
        # PCA
        fig, axs = plt.subplots(1,4,figsize=(16,4))

        scatter(df_, x, y, by='sample', c=colors['sample'],ax=axs[0], s=100)
        format_ax(ax=axs[0], xlabel=x, ylabel=y)
        add_legend(label='Sample', colors=colors['sample'], ax=axs[0], 
                loc='lower left', bbox_to_anchor=(0,0))
        scatter(df_, x, y, by='condition', c=colors['condition'], ax=axs[1], s=100)
        format_ax(ax=axs[1], xlabel=x, ylabel=y)
        add_legend(label='Condition', colors=colors['condition'], ax=axs[1], 
                loc='lower left', bbox_to_anchor=(0,0))
        scatter(df_, x, y, by='origin', c=colors['origin'], ax=axs[2], s=100)
        format_ax(ax=axs[2], xlabel=x, ylabel=y)
        add_legend(label='Origin', colors=colors['origin'], ax=axs[2], 
                loc='lower left', bbox_to_anchor=(0,0))
        scatter(df_, x, y, by='dataset', c=colors['dataset'], ax=axs[3], s=100)
        format_ax(ax=axs[3], xlabel=x, ylabel=y)
        add_legend(label='Dataset', colors=colors['dataset'], ax=axs[3], 
                loc='lower left', bbox_to_anchor=(0,0))

        fig.tight_layout()
        fig.savefig(
            os.path.join(path_results, f'{x}_{y}_pcs.png'),
            dpi=500
        )

# GSEA
for i in range(1,6):
    fig = PCA_gsea_loadings_plot(
        loads, adata.var, organism='human', collection='MSigDB_Hallmark_2020', i=i
    )
    fig.savefig(
        os.path.join(path_results, f'PC{i}_GSEA.png'),
        dpi=500
    )


##
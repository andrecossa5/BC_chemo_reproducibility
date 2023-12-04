"""
Correlation sc/bulk
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr
from plotting_utils._utils import *
from plotting_utils._colors import *
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
matplotlib.use('macOSX')


##


def reverse_complement(x):
    d = {'A':'T','C':'G','T':'A','G':'C'}
    L = list(x)
    new_L = []
    for x in L:
        new_L.append(d[x])
    rev = ''.join(new_L[::-1])
    return rev


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal')

# Read cells meta
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
meta = adata.obs[['GBC', 'sample', 'origin', 'condition', 'dataset']]

# sc freqs only PTs
df_freq = (
    meta.query('origin == "PT"')
    .groupby(['condition', 'origin', 'sample'])
    ['GBC'].value_counts(normalize=True).loc[lambda x:x>0]
    .reset_index(name='freq_sc')
    .rename(columns={'level_3':'GBC'})
    .merge(
        meta[['sample', 'dataset']]
        .reset_index(drop=True)
        .drop_duplicates(),
        
    )
    .drop_duplicates()
)

# Bulk freqs
df_bulk = pd.read_csv(
    os.path.join(path_data, 'perturb_bulk', 'bulk_prevalences.csv'), 
    index_col=0
)

# Prep freq df
df_ = (
    df_freq
    .assign(GBC_bulk=lambda x: x['GBC'].map(lambda x: reverse_complement(x)))
    .query('GBC_bulk in @df.index and freq_sc >= .001')
    .set_index(['GBC_bulk', 'sample'])[['freq_sc']]
    .join(
        df_bulk
        .reset_index().rename(columns={'index':'GBC_bulk'})
        .set_index(['GBC_bulk', 'sample'])[['cellular_prevalence_wi']]
        .rename(columns={'cellular_prevalence_wi':'freq_bulk'})
        .query('freq_bulk >= .001'),
        how='inner'
    )
)

# Corr plot
fig, ax = plt.subplots(figsize=(5,5))
scatter(df_.query('freq_bulk<.3 and freq_sc<.3'), 'freq_bulk', 'freq_sc', c='k', ax=ax, s=15)
sns.regplot(data=df_.query('freq_bulk<.3 and freq_sc<.3'), x='freq_bulk', y='freq_sc', ax=ax, scatter=False)
pho, p = pearsonr(df_["freq_bulk"], df_["freq_sc"])
format_ax(
    ax, title=f'Bulk-sc correlation: {pho:.2f} (p={p:.2f}, n={df_.shape[0]})',
    xlabel='Bulk prevalence', ylabel='Single-cell prevalence', reduce_spines=True
)
fig.tight_layout()

fig.savefig(os.path.join(path_results, 'bulk_sc_correlation.png'), dpi=300)


##
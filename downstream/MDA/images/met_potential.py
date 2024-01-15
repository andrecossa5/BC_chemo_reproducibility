"""
Clonal summary stats plots.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from plotting_utils._utils import *
from plotting_utils._colors import *
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
from BC_chemo_utils.clonal_utils import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal')

# Read cells meta
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
meta = meta[['GBC', 'sample', 'origin', 'condition', 'dataset']]

# sc freqs 
df_freq = (
    meta[['GBC', 'sample']]
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='freq')
    .rename(columns={'level_3':'GBC'})
    .merge(meta[['sample', 'origin']].drop_duplicates(), on='sample')    
)

# Expansions
df = (
    df_freq
    .merge(meta[['sample', 'dataset']].drop_duplicates(), on='sample')
    .pivot_table(index='GBC', columns=['dataset', 'origin'], values='freq')
    .melt(value_name='freq', ignore_index=False)
    .reset_index()
    .pivot_table(index=['GBC', 'dataset'], columns='origin', values='freq')
    .reset_index().dropna().set_index('GBC')
    .assign(
        met_potential=lambda x: x['lung']/x['PT'],
        logPT=lambda x: np.log2(x['PT']),
        loglung=lambda x: np.log2(x['lung'])
    )
)

fig, ax = plt.subplots(figsize=(5,4.5))
scatter(df, 'logPT', 'loglung', s=100, by='met_potential', c='Spectral_r', ax=ax)
format_ax(
    title=f'Clonal expansions, PT-lung (n={df.shape[0]})',
    xlabel='PT prevalence',
    ylabel='lung prevalence',
    reduce_spines=True,
    ax=ax
)
add_cbar(
    df['met_potential'], palette='Spectral_r', ax=ax, label='met_potential',
)
fig.tight_layout()
plt.show()













# fig.savefig(os.path.join(path_results, 'met_potential.png'), dpi=300)


##
















# n longitudinal with more than <x> cells
df = df.reset_index()


# Change...
L = []
origin = "lung"
for i in range(df.shape[0]):
    d = df.iloc[i,:].to_dict()
    dataset = d['dataset']
    GBC = d['GBC']
    L.append(
        df.query('origin==@origin and dataset==@dataset and GBC==@GBC').shape[0]
    )
df['n_lung'] = L

# Longitudinal clones
(
    df
    .query('n_PT>=10 and n_lung>=10')
    .sort_values('met_potential', ascending=False)
    .to_csv(os.path.join(path_data, 'longitudinal_clones.csv'), index=False)
)


##
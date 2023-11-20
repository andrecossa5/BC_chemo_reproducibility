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
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
meta = adata.obs[['GBC', 'sample', 'origin', 'condition', 'dataset']]

# sc freqs 
df_freq = (
    meta.groupby(['condition', 'origin', 'sample'])
    ['GBC'].value_counts(normalize=True).loc[lambda x:x>0]
    .reset_index(name='freq')
    .rename(columns={'level_3':'GBC'})
    .merge(
        meta[['sample', 'dataset']]
        .reset_index(drop=True)
        .drop_duplicates(),
        
    )
    .drop_duplicates()
)


##


# Cum clone percentages
fig, axs = plt.subplots(1,2,figsize=(9.5,5))

for sample in df_freq.query('origin == "PT"')['sample'].unique():
    color = 'grey' if sample.startswith('AC') else 'k'
    print(color)
    x = df_freq.query('sample == @sample')['freq'].cumsum().values
    axs[0].plot(x, '.-', color=color)
format_ax(axs[0], title='PT (n=12)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

for sample in df_freq.query('origin == "lung"')['sample'].unique():
    color = 'grey' if sample.startswith('NT_NT') else 'k'
    print(color)
    x = df_freq.query('sample == @sample')['freq'].cumsum().values
    axs[1].plot(x, '.-', color=color)
format_ax(axs[1], title='lung (n=12)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

add_legend(
    label='Sample type', colors={'treated':'grey', 'untreated':'k'}, ax=axs[0],
    bbox_to_anchor=(.95,.05), loc='lower right', ticks_size=8, label_size=10, artists_size=8
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cum_fractions.png'), dpi=300)


##


# Prep stats 
df_stats = stats_summary(df_freq, freq='freq')
df_stats = (
    df_freq[['sample', 'origin', 'dataset', 'condition']]
    .merge(df_stats, on='sample')
    .drop_duplicates()
)
df_stats['condition'] = df_stats['condition'].astype('str')
df_stats['origin'] = df_stats['origin'].astype('str')

# Reformat condition
tests = [
    (df_stats['condition'].isin(['AC_AC', 'AC_NT'])) & (df_stats['origin'] == 'PT'),
    (df_stats['condition'].isin(['NT_NT', 'NT_AC'])) & (df_stats['origin'] == 'PT'),
    (df_stats['condition'] == 'NT_NT') & (df_stats['origin'] == 'lung'),
    (df_stats['condition'] == 'AC_NT') & (df_stats['origin'] == 'lung'),
    (df_stats['condition'] == 'NT_AC') & (df_stats['origin'] == 'lung'),
    (df_stats['condition'] == 'AC_AC') & (df_stats['origin'] == 'lung')
]
labels = [
    'PT, treated', 'PT, untreated', 'lung, untreated', 
    'lung, PT-treated', 'lung, lung-treated',
    'lung, double-treated'
]
df_stats['condition'] = np.select(tests, labels)
df_stats['condition'] = pd.Categorical(
    df_stats['condition'], 
    categories=[
        'PT, treated', 'PT, untreated', 'lung, untreated',
        'lung, PT-treated', 'lung, lung-treated', 'lung, double-treated'
    ]
)


##


# Viz, boxplots
pairs = [
    ['PT, untreated', 'PT, treated'],
    ['lung, untreated', 'lung, PT-treated'],
    ['lung, untreated', 'lung, lung-treated'],
    ['lung, lung-treated', 'lung, double-treated'],
    ['lung, PT-treated', 'lung, double-treated'],
]

fig, axs = plt.subplots(2,1, figsize=(11,6), sharex=True)

box(df_stats, 'condition', 'n', ax=axs[0], c='white', with_stats=True, pairs=pairs)
strip(df_stats, 'condition', 'n', ax=axs[0], c='k')
format_ax(ax=axs[0], title='n clones by condition and origin', 
        ylabel='n clones', reduce_spines=True)

box(df_stats, 'condition', 'SH', ax=axs[1], c='white', with_stats=True, pairs=pairs)
strip(df_stats, 'condition', 'SH', ax=axs[1], c='k')
format_ax(ax=axs[1], title='Shannon Entropy (SH) by condition and origin', ylabel='SH', reduce_spines=True)

fig.subplots_adjust(hspace=.4)
fig.savefig(os.path.join(path_results, 'stats.png'), dpi=300)


##


# Viz expansions
df_ = (
    df_freq
    .pivot_table(index='GBC', columns=['dataset', 'origin'], values='freq')
    .melt(value_name='freq', ignore_index=False)
    .reset_index()
    .pivot_table(index=['GBC', 'dataset'], columns='origin', values='freq')
    .reset_index().dropna().set_index('GBC')
    .assign(met_potential=lambda x: x['lung']/x['PT'])
)

fig, ax = plt.subplots(figsize=(5,4.5))
scatter(df_, 'PT', 'lung', s=100, by='met_potential', c='Spectral_r', ax=ax)
format_ax(
    title=f'Clonal expansions, PT-lung (n={df_.shape[0]})',
    xlabel='PT prevalence',
    ylabel='lung prevalence',
    reduce_spines=True,
    ax=ax
)
add_cbar(
    df_['met_potential'], palette='Spectral_r', ax=ax, label='met_potential',
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'met_potential.png'), dpi=300)


##
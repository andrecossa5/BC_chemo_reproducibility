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
df = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
df = df[['GBC', 'sample', 'origin', 'condition', 'dataset']]

# sc freqs 
df_freq = (
    df[['GBC', 'sample']]
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='freq')
    .rename(columns={'level_3':'GBC'})
    .merge(df[['sample', 'origin']].drop_duplicates(), on='sample')    
)


##


# Cum clone percentages
fig, axs = plt.subplots(1,2,figsize=(9.5,5))

for sample in df_freq.query('origin == "PT"')['sample'].unique():
    color = 'grey' if sample.startswith('AC') else '#2A6CB6'
    x = df_freq.query('sample == @sample')['freq'].cumsum().values
    axs[0].plot(x, '-', color=color, linewidth=2)
format_ax(axs[0], title='PT (n=12)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

for sample in df_freq.query('origin == "lung"')['sample'].unique():
    color = 'grey' if sample.startswith('NT_NT') else '#2A6CB6'
    x = df_freq.query('sample == @sample')['freq'].cumsum().values
    axs[1].plot(x, '-', color=color, linewidth=2)
format_ax(axs[1], title='lung (n=9)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

add_legend(
    label='Sample type', colors={'treated':'#2A6CB6', 'untreated':'grey'}, ax=axs[0],
    bbox_to_anchor=(.95,.05), loc='lower right', ticks_size=8, label_size=10, artists_size=8
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cum_fractions.png'), dpi=300)


##


# Prep stats 
df_stats = stats_summary(df_freq, freq='freq')
df_stats = (
    df_freq[['sample', 'origin']]
    .merge(df[['sample', 'condition']], on='sample').drop_duplicates()
    .merge(df_stats, on='sample').drop_duplicates()
)
df_stats['condition'] = df_stats['condition'].astype('str')
df_stats['origin'] = df_stats['origin'].astype('str')

# Re-set categories
df_stats['condition'] = pd.Categorical(
    df_stats['condition'], 
    categories=[
        'PT, treated', 'PT, untreated', 'lung, untreated',
        'lung, single-treated', 'lung, double-treated'
    ]
)


##


# Viz, boxplots
pairs = [
    ['PT, untreated', 'PT, treated'],
    ['lung, untreated', 'lung, single-treated'],
    ['lung, single-treated', 'lung, double-treated'],
    ['lung, untreated', 'lung, double-treated'],
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
    .merge(df[['sample', 'dataset']].drop_duplicates(), on='sample')
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


# Longitudinal with more than 10 cells
df_ = (
    df.groupby(["GBC", "dataset", "origin"]).size().reset_index(name="cell_count")
    .assign(frequency=lambda x: 
        x.groupby(["dataset", "origin"])["cell_count"].transform(lambda x: x / x.sum())
    )
)
df_ = df_.pivot_table(columns=['origin'], index=['GBC', 'dataset'], values='cell_count', fill_value=0).reset_index().merge(
    df_.pivot_table(columns=['origin'], index=['GBC', 'dataset'], values='frequency', fill_value=0).reset_index(),
    on=['GBC', 'dataset'], suffixes=['_count', '_freq']
)
df_ = df_.query('PT_count>=10 and lung_count>=10')
df_.to_csv(os.path.join(path_data, 'longitudinal_clones.csv'))


##

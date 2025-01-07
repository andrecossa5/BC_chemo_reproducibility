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
path_data = os.path.join(path_main, 'data', 'CTCs')
path_results = os.path.join(path_main, 'results', 'CTCs', 'clonal')

# Read cells meta
df = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
df = df[['GBC', 'sample', 'origin', 'dataset']]

# sc freqs 
df_freq = (
    df[['GBC', 'sample']].groupby('sample')['GBC'].value_counts(normalize=True).to_frame('freq')
    .join(
        df[['GBC', 'sample']].groupby('sample')['GBC'].value_counts().to_frame('n')
    )
    .reset_index()
    .merge(df[['sample', 'origin']].drop_duplicates(), on='sample')    
)
df_freq['origin'] = pd.Categorical(df_freq['origin'], categories=['PT', 'CTC', 'lung'])


##


# Cum clone percentages
fig, ax = plt.subplots(figsize=(5,5))

colors = create_palette(df_freq, 'origin', ten_godisnot)

for sample in df_freq['sample'].unique():
    df_ = df_freq.query('sample == @sample')
    origin = df_['origin'].unique()[0]
    x = df_freq.query('sample == @sample')['freq'].cumsum().values
    ax.plot(x, '-', color=colors[origin], linewidth=2)
    format_ax(ax, title='Cumulative clone percentages (n=12 samples)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

add_legend(
    label='Sample type', colors=colors, ax=ax,
    bbox_to_anchor=(.95,.05), loc='lower right', ticks_size=8, label_size=10, artists_size=8
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cum_fractions.png'), dpi=300)


##


# Prep stats 
df_stats = stats_summary(df_freq, freq='freq')
df_stats = (
    df_freq[['sample', 'origin']]
    .merge(df[['sample', 'dataset']], on='sample').drop_duplicates()
    .merge(df_stats, on='sample').drop_duplicates()
)
df_stats.to_csv(os.path.join(path_results, 'summary_clonal_stats.csv'))

# Re-set categories
df_stats['origin'] = pd.Categorical(df_stats['origin'], categories=['PT', 'CTC', 'lung'])
df_stats['dataset'] = pd.Categorical(df_stats['dataset'], categories=['1', '2', '3', '4'])


##


# Viz, boxplots
pairs = [
    ['PT', 'CTC'],
    ['PT', 'lung'],
    ['lung', 'CTC']
]

fig, axs = plt.subplots(1,2,figsize=(8,4))

df_stats['n'] = np.log10(df_stats['n'])

ax = axs[0]
box(df_stats, 'origin', 'n', ax=ax, c='white', with_stats=True, pairs=pairs)
strip(df_stats, 'origin', 'n', ax=ax, c='k')
format_ax(ax=axs[0], title='n clones by origin', 
        ylabel='n clones', reduce_spines=True)

ax = axs[1]
box(df_stats, 'origin', 'SH', ax=ax, c='white', with_stats=True, pairs=pairs)
strip(df_stats, 'origin', 'SH', ax=ax, c='k')
format_ax(ax=axs[1], title='Shannon Entropy (SH) by origin', ylabel='SH', reduce_spines=True)

fig.subplots_adjust(hspace=.4)
fig.savefig(os.path.join(path_results, 'stats.png'), dpi=300)


##


# n longitudinal
df_long = df_freq.merge(df[['sample', 'dataset']], on='sample', how='left').drop_duplicates()
freq = df_long.pivot_table(index=['dataset', 'GBC'], columns='origin', values='freq').rename(columns={'PT':'freq_PT', 'CTC':'freq_CTC', 'lung':'freq_lung'})
n = df_long.pivot_table(index=['dataset', 'GBC'], columns='origin', values='n').rename(columns={'PT':'n_PT', 'CTC':'n_CTC', 'lung':'n_lung'})
df_long = freq.join(n)
df_long = df_long.dropna().reset_index()

# df_long
# np.sum(np.sum(df_long.iloc[:,-3:]>=10, axis=1)==3)

# Save
df_long.to_csv(os.path.join(path_results, 'logitudinal_clones.csv'))




##


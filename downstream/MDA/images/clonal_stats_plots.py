"""
Clonal summary stats plots.
"""

import os
import numpy as np
import pandas as pd
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
meta = meta[['GBC', 'sample', 'origin', 'condition', 'dataset']].rename(columns={'condition':'condition_sample'})

# Remove AC_NT branch
meta = meta.loc[~meta['dataset'].str.contains('AC_NT'),:].copy()

# Calculate clone cell count and frequency
df_freq = (
    meta
    .groupby(["GBC", "dataset", "origin"])
    .size().reset_index(name="cell_count")
    .assign(
        freq=lambda x: x.groupby(["dataset", "origin"])["cell_count"].transform(lambda x: x / x.sum())
    )
    .merge(meta[['sample', 'dataset', 'origin', 'condition_sample']].drop_duplicates(), on=['dataset', 'origin'])  
    [['GBC', 'sample', 'dataset', 'origin', 'condition_sample', 'cell_count', 'freq']]
)


##


#1. Cum clone percentages
fig, axs = plt.subplots(1,2,figsize=(8,4))

for sample in df_freq.query('origin == "PT"')['sample'].unique():
    color = 'grey' if sample.startswith('NT') else '#2A6CB6'
    x = df_freq.query('sample == @sample')['freq'].sort_values(ascending=False).cumsum().values
    axs[0].plot(x, '-', color=color, linewidth=2)
format_ax(axs[0], title='PT (n=9)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

for sample in df_freq.query('origin == "lung"')['sample'].unique():
    color = 'grey' if sample.startswith('NT') else '#2A6CB6'
    x = df_freq.query('sample == @sample')['freq'].sort_values(ascending=False).cumsum().values
    axs[1].plot(x, '-', color=color, linewidth=2)
format_ax(axs[1], title='lung (n=9)', xlabel='Clones, ranked', ylabel='Cumulative fraction')

add_legend(
    label='Sample type', colors={'treated':'#2A6CB6', 'untreated':'grey'}, ax=axs[0],
    bbox_to_anchor=(.95,.05), loc='lower right', ticks_size=8, label_size=10, artists_size=8
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cum_fractions.png'), dpi=300)


##


#2. Per-sample clonal stats, grouped by condition_sample
df_stats_sample = stats_summary(df_freq, freq='freq')
df_stats_sample = df_freq.iloc[:,1:-2].drop_duplicates().merge(df_stats_sample.reset_index(), on='sample')
df_stats_sample['condition_sample'] = pd.Categorical(
    df_stats_sample['condition_sample'], 
    categories=[
        'PT, treated', 'PT, untreated', 'lung, untreated',
        'lung, single-treated', 'lung, double-treated'
    ]
)
pairs = [
    ['PT, untreated', 'PT, treated'],
    ['lung, untreated', 'lung, single-treated'],
    ['lung, single-treated', 'lung, double-treated'],
    ['lung, untreated', 'lung, double-treated'],
]

fig, axs = plt.subplots(2,1, figsize=(10,6), sharex=True)

box(df_stats_sample, 'condition_sample', 'n', ax=axs[0], c='white', with_stats=True, pairs=pairs)
strip(df_stats_sample, 'condition_sample', 'n', ax=axs[0], c='k')
format_ax(ax=axs[0], title='n clones by condition and origin', ylabel='n clones', reduce_spines=True)

box(df_stats_sample, 'condition_sample', 'SH', ax=axs[1], c='white', with_stats=True, pairs=pairs)
strip(df_stats_sample, 'condition_sample', 'SH', ax=axs[1], c='k')
format_ax(ax=axs[1], title='Shannon Entropy (SH) by condition and origin', ylabel='SH', reduce_spines=True)

fig.subplots_adjust(hspace=.4)
fig.savefig(os.path.join(path_results, 'samples_stats.png'), dpi=300)


##


#3. Metastatic potential
df_long = pd.merge(
    df_freq.pivot_table(index=['GBC', 'dataset'], columns=['origin'], values='cell_count', fill_value=0).reset_index(),
    df_freq.pivot_table(index=['GBC', 'dataset'], columns=['origin'], values='freq', fill_value=0).reset_index()
    , on=['GBC', 'dataset'], suffixes=['_count', '_freq'] 
)
df_long['met_potential'] = df_long['lung_freq'] / (df_long['PT_freq'])

fig, ax = plt.subplots(figsize=(4.5,4.5))
scatter(df_long.query('PT_freq>0 and lung_freq>0'), 'PT_freq', 'lung_freq', s=100, by='met_potential', c='Spectral_r', ax=ax)
format_ax(
    title=f'Clonal expansions, PT-lung (n={df_long.shape[0]})',
    xlabel='PT prevalence',
    ylabel='lung prevalence',
    reduce_spines=True,
    ax=ax
)
add_cbar(
    df_long['met_potential'], palette='Spectral_r', ax=ax, label='met_potential',
)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'met_potential.png'), dpi=300)


##


#2. Per-dataset clonal stats, grouped by condition_dataset
df_long['condition_dataset'] = np.select([
    df_long['dataset'].str.startswith('NT_NT'),
    df_long['dataset'].str.startswith('NT_AC'),
    df_long['dataset'].str.startswith('AC_AC')
    ], 
    ['Non treated', 'Adj', 'Neo + Adj']
)


##


# Per mice number
x = 1

# n total clones in PT (n cells>=min n cells)
df_summary = (
    df_long.query('PT_count>=@x')
    .groupby(['condition_dataset', 'dataset'])['GBC'].nunique()
    .to_frame('n total clones (PT)')
)

# n total clones in lung (n cells>=min n cells)
df_summary = df_summary.join(
    df_long.query('lung_count>=@x')
    .groupby(['condition_dataset', 'dataset'])['GBC'].nunique()
    .to_frame('n total clones (lung)')
)

# n total clones (n cells>=min n cells)
df_summary = df_summary.join(
    df_long.query('PT_count>=@x or lung_count>=@x')
    .groupby(['condition_dataset', 'dataset'])['GBC'].nunique()
    .to_frame('n total clones')
)

# n prometastatic clones (n cells>=min n cells, both organs)
df_summary = df_summary.join(
    df_long.query('PT_count>=@x and lung_count>=@x')
    .groupby(['condition_dataset', 'dataset'])['GBC'].nunique()
    .to_frame('n prometastatic clones')
)

# % prometastatic clones in PT (n prometastatic / n total, as defined above)
df_summary['% prometastatic clones (PT)'] = df_summary['n prometastatic clones'] / df_summary['n total clones (PT)']
# df_summary['% prometastatic clones, min-max'] = (df_summary['% prometastatic clones'] - df_summary['% prometastatic clones'].min()) / \
#                                                 ( df_summary['% prometastatic clones'].max() - df_summary['% prometastatic clones'].min() )

# % prometastatic cells in PT (n cells from prometastatic clones / n total cells, as defined above)
df_summary['% prometastatic cells (PT)'] = df_long.query('PT_count>=@x and lung_count>=@x').groupby(['condition_dataset', 'dataset'])['PT_count'].sum() / \
                                      df_long.query('PT_count>=1').groupby(['condition_dataset', 'dataset'])['PT_count'].sum()
# df_summary['% prometastatic cells, min-max'] = (df_summary['% prometastatic cells'] - df_summary['% prometastatic cells'].min()) / \
#                                                 ( df_summary['% prometastatic cells'].max() - df_summary['% prometastatic cells'].min() )

# Additional PTs clonal stats
df_PTs = df_stats_sample.query('origin=="PT"')[['dataset', 'min_freq', 'max_freq', 'median_freq', 'std_freq', 'SH']]
df_PTs.columns = ['dataset'] + [ f'{x} (PT)' for x in df_PTs.columns[1:] ]
df_summary = df_summary.reset_index().merge(df_PTs, on='dataset')

# Additional lungs clonal stats
df_lungs = df_stats_sample.query('origin=="lung"')[['dataset', 'min_freq', 'max_freq', 'median_freq', 'std_freq', 'SH']]
df_lungs.columns = ['dataset'] + [ f'{x} (lung)' for x in df_lungs.columns[1:] ]
df_summary = df_summary.merge(df_lungs, on='dataset')

df_summary = df_summary.set_index(['condition_dataset', 'dataset']).T
df_summary.to_csv(os.path.join(path_results, 'clonal_stats_summary.csv'))
df_summary


##


# Other 
t = [0,1,2,5,10]


# n prometastatic clones
d = {}
for x in t:
    d[x] = df_long.query('PT_count>@x and lung_count>@x').groupby('condition_dataset')['GBC'].nunique().to_dict()
df_plot = pd.DataFrame(d).reset_index(names='condition_dataset').melt(id_vars=['condition_dataset'], var_name='min n cell', value_name='n pro-metastatic clones')
df_plot['condition_dataset'] = pd.Categorical(df_plot['condition_dataset'], categories=['Non treated', 'Adj', 'Neo + Adj'])

fig, ax = plt.subplots(figsize=(4,4))
colors = create_palette(df_plot, 'condition_dataset', 'Reds') 
sns.barplot(data=df_plot, x='min n cell', y='n pro-metastatic clones', hue='condition_dataset', palette=colors, ax=ax)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, 'n pro-metastatic clones.png'), dpi=300)

# % of prometastatic clones
d = {}
for x in t:
    n = df_long.query('PT_count>@x and lung_count>@x').groupby('condition_dataset')['GBC'].nunique() / \
        df_long.query('PT_count>0').groupby('condition_dataset')['GBC'].nunique()
    d[x] = np.round(n,2).to_dict()
df_plot = pd.DataFrame(d).reset_index(names='condition_dataset').melt(id_vars=['condition_dataset'], var_name='min n cell', value_name='% pro-metastatic clones')
df_plot['condition_dataset'] = pd.Categorical(df_plot['condition_dataset'], categories=['Non treated', 'Adj', 'Neo + Adj'])

fig, ax = plt.subplots(figsize=(4,4))
colors = create_palette(df_plot, 'condition_dataset', 'Reds') 
sns.barplot(data=df_plot, x='min n cell', y='% pro-metastatic clones', hue='condition_dataset', palette=colors, ax=ax)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, '% pro-metastatic clones.png'), dpi=300)
 
# % of prometastatic clones, aggregated by mice and condition
L = []
for x in t:
    n = df_long.query('PT_count>@x and lung_count>@x').groupby(['condition', 'dataset'])['GBC'].nunique() / \
        df_long.query('PT_count>0').groupby(['condition', 'dataset'])['GBC'].nunique()
    n = n.loc[lambda x: ~x.isna()]
    n = n.reset_index().assign(min_n_cell=x).rename(columns={'GBC':'% pro-metastatic clones', 'min_n_cell':'min n cell'})
    L.append(n)
df_plot = pd.concat(L)
df_plot['condition'] = pd.Categorical(df_plot['condition'], categories=['Non treated', 'Adjuvant treated', 'Double treated'])

df_plot.loc[df_plot['min n cell']==10]

fig, ax = plt.subplots(figsize=(4,4))
colors = create_palette(df_plot, 'condition', 'Reds') 
sns.barplot(data=df_plot, x='min n cell', y='% pro-metastatic clones', hue='condition', palette=colors, ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_results, '% pro-metastatic clone by mice.png'), dpi=300)

# # % of cells from prometastatic clones in the PT, aggregated by mice and condition
# L = []
# for x in t:
#     n = df_.query('PT_count>@x and lung_count>@x').groupby(['condition', 'dataset'])['PT_count'].sum() / \
#         df_.query('PT_count>0').groupby(['condition', 'dataset'])['PT_count'].sum()
#     n = n.loc[lambda x: ~x.isna()]
#     n = n.reset_index().assign(min_n_cell=x).rename(columns={'PT_count':'% pro-metastatic cells', 'min_n_cell':'min n cell'})
#     L.append(n)
# df_plot = pd.concat(L)
# df_plot['condition'] = pd.Categorical(df_plot['condition'], categories=['Non treated', 'Adjuvant treated', 'Double treated'])
# 
# df_plot.loc[df_plot['min n cell']==10]
# 
# fig, ax = plt.subplots(figsize=(4,4))
# colors = create_palette(df_plot, 'condition', 'Reds') 
# sns.barplot(data=df_plot, x='min n cell', y='% pro-metastatic cells', hue='condition', palette=colors, ax=ax)
# fig.tight_layout()
# fig.savefig(os.path.join(path_results, '% pro-metastatic cells by mice.png'), dpi=300)
# 
# 
# ##
# 

# df_ = df_.query('PT_count>=10 and lung_count>=10')
# df_.to_csv(os.path.join(path_data, 'longitudinal_clones.csv'))


##

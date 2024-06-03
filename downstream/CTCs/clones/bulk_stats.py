"""
Essential bulk analysis script.
"""

import os
import numpy as np
import pandas as pd
from plotting_utils._plotting import *


## 


# Paths
run = 'bulk_treatment'
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'CTCs', run)      # Or other bulk output folder
path_results = os.path.join(path_main, 'results', 'CTCs', 'clonal')

# Read prevalences
df = pd.read_csv(os.path.join(path_data, 'bulk_GBC_reference.csv'), index_col=0)
common = pd.read_csv(os.path.join(path_data, 'common.csv'), index_col=0)

# Reformat
tests = [ df['sample'].str.contains('PT'), df['sample'].str.contains('LN'), df['sample'].str.contains('CTC') ] 
df['origin'] = np.select(tests, ['PT', 'LN', 'CTC'], default='Ref')
df['dataset'] = df['sample'].map(lambda x: '_'.join(x.split('_')[-1:]))
tests = [ df['sample'].str.contains('NT'), df['sample'].str.contains('AC'), df['sample'].str.contains('line') ] 
df['treatment'] = np.select(tests, ['NT', 'AC', 'line'], default=np.nan)
df['var'] = df['treatment'] + '_' + df['origin'] 


##


# n clones per sample
df_ = (
    df.groupby('sample')
    .apply(lambda x: x.index.unique().size)
    .sort_values(ascending=False)
    .to_frame('n')
    .reset_index()
)
fig, ax = plt.subplots(figsize=(8,4.5))
bar(df_, 'n', 'sample', s=.75, c='k', a=.7, ax=ax)
format_ax(ax=ax, title='n clones by sample', ylabel='n', xticks=df_['sample'], rotx=90)
ax.spines[['left', 'top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'n_clones_{run}.png'), dpi=500)

# Per origin
df_ = df_.merge(df[['sample', 'var']], on='sample').drop_duplicates().set_index('sample')

fig, ax = plt.subplots(figsize=(4,4))
order = ['NT_PT', 'AC_PT', 'NT_CTC', 'AC_CTC', 'NT_LN']
box(
    df_, 
    x='var', y='n', ax=ax, c='grey', with_stats=True, 
    pairs=[['NT_PT', 'AC_PT'], ['AC_CTC', 'NT_CTC']], 
    order=order
)
strip(df_, x='var', y='n', ax=ax, c='k', order=order)
format_ax(ax=ax, title='n clones', ylabel='n', rotx=90, reduce_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'n_clones_by_origin_and_treatment_{run}.png'), dpi=300)


##


# Cumulative clone percentage, all samples
colors = create_palette(df, 'var', ten_godisnot)

fig, ax = plt.subplots(figsize=(4.5,4.5))
for s in df['sample'].unique():
    df_ = df.query('sample==@s')
    x = (df_['read_count'] / df_['read_count'].sum()).cumsum()
    var = df.query('sample==@s')['var'].unique()[0]
    ax.plot(range(len(x)), x, c=colors[var], linewidth=2.5)

ax.set(title='Clone prevalences', xlabel='Ranked clones', ylabel='Cumulative frequence')
add_legend(ax=ax, colors=colors, bbox_to_anchor=(1,0), loc='lower right', ticks_size=8, label_size=10, artists_size=8)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'cum_percentages_{run}.png'), dpi=300)


##


# Shannon entropies
SH = []
for s in df['sample'].unique():
    df_ = df.query('sample==@s')
    x = df_['read_count'] / df_['read_count'].sum()
    SH.append(-np.sum( np.log10(x) * x ))
df_ = (pd.Series(SH, index=df['sample'].unique())
    .to_frame('SH')
    .sort_values(by='SH', ascending=False)
    .reset_index().rename(columns={'index':'sample'})
    .merge(df[['sample', 'var']], on='sample')
    .drop_duplicates()
    .set_index('sample')
)
fig, ax = plt.subplots(figsize=(4,4))
order = ['NT_PT', 'AC_PT', 'NT_CTC', 'AC_CTC', 'NT_LN']
box(
    df_, 
    x='var', y='SH', ax=ax, c='grey', with_stats=True, 
    pairs=[['NT_PT', 'AC_PT'], ['AC_CTC', 'NT_CTC']], 
    order=order
)
strip(df_, x='var', y='SH', ax=ax, c='k', order=order)
format_ax(ax=ax, title='Shannon entropy', ylabel='SH', rotx=90, reduce_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'SH_by_origin_and_treatment_{run}.png'), dpi=300)


##


# Commmon clones
fig, ax = plt.subplots(figsize=(5.5,5))
ax = sns.heatmap(data=common, ax=ax, robust=True, cmap='mako', annot=True, xticklabels=True, 
    yticklabels=True, annot_kws={'size':5}, cbar=True, fmt='.0f',
    cbar_kws={'fraction':0.05, 'aspect':35, 'pad': 0.02, 'shrink':1.0, 'label':'n clones'}
)
ax.set(title='n common clones', xlabel='Sample', ylabel='Sample')
ax.set_xticklabels(ax.get_xticklabels(), fontsize=7)
ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'common_{run}.png'), dpi=300)


##


# Fish, longitudinal
# ...

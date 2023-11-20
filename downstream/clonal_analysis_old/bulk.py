"""
Bulk clonal analysis with pipeline nf-pgp-perturbseq 23-10-2023
"""

import os
import pandas as pd
import numpy as np
from plotting_utils._plotting_base import *
from Cellula._utils import rescale
from BC_chemo_utils.clonal_utils import *
import matplotlib
matplotlib.use('MacOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'MDA', 'perturb_bulk')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal', 'bulk')

# Read data
# L = []
# for root_path, _, files in os.walk(path_data):
#     for f in files:
#         if f == 'clonal_prevalences.csv':
#             sample = root_path.split('/')[-1]
#             L.append(pd.read_csv(
#                 os.path.join(root_path, 'clonal_prevalences.csv'), index_col=0)
#                 .assign(sample=sample)
#             )
# df = pd.concat(L)
# df.to_csv(os.path.join(path_data, 'bulk_prevalences.csv'))

df = pd.read_csv(os.path.join(path_data, 'bulk_prevalences.csv'), index_col=0)

# Filter only clones found with both approaches
df = df.loc[df['found_wo'] & df['found_wi']].copy()

# Colors
condition_colors = create_palette(df, 'condition', ['g', 'y'])


##


# Fig n clones
#======================================================================#

fig, axs = plt.subplots(1, 2, figsize=(8, 5), sharey=True)

df_ = df.groupby('sample').size().to_frame(name='n')
bar(df_, 'n', c='k', s=0.8, a=0.8, ax=axs[0], annot_size=8)
format_ax(axs[0], title='n clones by sample', ylabel='n clones', 
    xticks=df_.index, rotx=90, log=True, reduce_spines=True
)

df_ = df.groupby('condition').size().to_frame(name='n').reset_index()
bar(df_, 'n', c='k', s=0.8, a=0.8, ax=axs[1], annot_size=10)
format_ax(axs[1], title='n clones by condition', 
    xticks=df_['condition'].unique(), log=True, reduce_spines=True
)

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'n_clones.png'), dpi=300)

#======================================================================#


##


# Summaries
df['sample'] = pd.Categorical(df['sample'])
df['condition'] = pd.Categorical(df['condition'])
sample_df = stats_summary(df, 'sample')
sample_df.to_csv(os.path.join(path_results, 'summary_sample.csv'))
condition_df = stats_summary(df, 'condition')
sample_df.to_csv(os.path.join(path_results, 'summary_condition.csv'))

d, counts, common_idxs = common(df) # 0.0, 3.71 median/mean
d, counts, common_idxs = common(df.query('condition == "untreated"')) # 0 5.37 median/mean
d, counts, common_idxs = common(df.query('condition == "treated"')) # 0 3.71 median/mean


##


# Fig n_commons
#======================================================================#

# All samples
counts = { 
    n : 
    common_in(df, n, normalize=False) \
    for n in range(2, df['sample'].unique().size) # n samples
}

df_ = pd.DataFrame(data=counts.values(), index=counts.keys(), columns=['n'])

fig, ax = plt.subplots(figsize=(6, 5))
bar(df_, 'n', c='k', s=0.8, a=0.8, ax=ax, annot_size=10)
format_ax(ax, title='n common clones', ylabel='n common clones', 
    xlabel='n tumors', xticks=df_.index
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'n_commons.png'), dpi=300)

#======================================================================#


##


# Fig difference in common clones distribution
#======================================================================#

# All samples
d, counts, common_idxs = common(df)

# Tests
combos = [ x.split('|') for x in common_idxs.keys() ]
counts_with_labels = { k : len(common_idxs[k]) for k in common_idxs }
df_ = (
    pd.DataFrame(
        data=counts_with_labels.values(), 
        columns=['n'])
    ).assign(
    between=counts_with_labels.keys()
)
df_[['a', 'b']] = df_['between'].str.split('|', expand=True)
test_tr = df_['a'].str.startswith('AC') & df_['b'].str.startswith('AC') 
test_un = df_['a'].str.startswith('NT') & df_['b'].str.startswith('NT') 
df_['pair_type'] = np.select([test_tr, test_un], ['treated', 'untreated'], default='mixed')
df_['pair_type'].value_counts()
df_['pair_type'] = pd.Categorical(df_['pair_type'], categories=['mixed', 'untreated', 'treated'])

# Fig
fig, ax = plt.subplots(figsize=(4.5, 5))
strip(df_, 'pair_type', 'n', c='k', ax=ax, with_stats=True, 
    pairs=[['treated', 'untreated'], ['treated', 'mixed']], 
    order=['mixed', 'untreated', 'treated']
)
format_ax(ax=ax, title='n common clones among samples', ylabel='n common', reduce_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'n_commons_per_condition.png'), dpi=300)

#======================================================================#


##

# Fig bubble
#======================================================================#

import random
n_clones = df.index.unique().size

clones_colors = { 
    clone : color for clone, color in \
    zip(
        df.index.unique(), 
        list(
            ''.join( ['#'] + [random.choice('ABCDEF0123456789') for i in range(6)] )  \
            for _ in range(n_clones)
        )
    )
}

df_ = df.reset_index().rename(columns={'index':'GBC'}).sample(df.shape[0])
df_['area_plot'] = rescale(df_['cellular_prevalence_wo']) * (3000-5) + 5

# Fig 
fig, ax = plt.subplots(figsize=(7, 7))
scatter(df_, 'GBC', 'sample', by='GBC', c=clones_colors, s='area_plot', a=0.5, ax=ax)
format_ax(ax, title='Clones by sample', xlabel='Clones', xticks='')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'bubbles.png'), dpi=300)

#======================================================================#


##


# Per sample metrics 
#======================================================================#

df_ = sample_df.copy()
df_['condition'] = df_.index.map(lambda x: 'treated' if x.startswith('AC') else 'untreated')
df_['condition'] = pd.Categorical(df_['condition'], categories=('untreated', 'treated'))
pairs = [('untreated', 'treated')]

# Fig
fig, axs = plt.subplots(1, 3, figsize=(7, 4))
box(df_, 'condition', 'std_freq', c='white', ax=axs[0], with_stats=True, pairs=pairs)
box(df_, 'condition', 'median_freq', c='white', ax=axs[1], with_stats=True, pairs=pairs)
box(df_, 'condition', 'SH', c='white', ax=axs[2], with_stats=True, pairs=pairs)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'metrics.png'), dpi=300)

plt.show()

#======================================================================#


##

# Distribution by condition
#======================================================================#

fig, ax = plt.subplots(figsize=(5, 5))

df_['Cell_fraction'] = -np.log10(df_['cellular_prevalence_wo'])

sns.kdeplot(data=df_, x='Cell_fraction', hue='condition', fill=True, ax=ax)
format_ax(ax, title='Clones frequency distribution', 
          xlabel='-log10(Cellular prevalence)', ylabel='Density')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'distributions.png'))

#======================================================================#


##


########################################################################


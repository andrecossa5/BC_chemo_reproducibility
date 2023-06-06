# Bulk

########################################################################

# Utilities
import sys
import os
import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
from statannotations.Annotator import Annotator
import seaborn as sns
import matplotlib.cm as cm

from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula._utils import *
matplotlib.use('MacOSX')

########################################################################

# Utils 


def format_df(df):
    '''Helper.'''
    cond_1 = np.array([ bool(re.search('C1|C3', x)) for x in df['sample'] ])
    cond_2 = np.array([ bool(re.search('C2|C4', x)) for x in df['sample'] ])
    cond_3 = df['sample'] == 'Reference'

    df['condition'] = np.select(
        [cond_1, cond_2, cond_3], 
        ['treated', 'untreated', 'in_vitro']
    )

    df['condition'] = pd.Categorical(
        df['condition'], categories=['in_vitro', 'untreated', 'treated']
    )

    df['sample'] = pd.Categorical(df['sample'])

    return df


##


def stats_summary(df, key='sample'):
    '''Basic stats summary.'''

    stats = {
        key : [],
        'n' : [],
        'min_freq' : [],
        'max_freq' : [],
        'median_freq' : [],
        'mean_freq' : [],
        'std_freq' : [],
        'SH' : [],
        'EI' : []
    }
    
    for c in df['sample'].cat.categories:

        df_one = df.query('sample == @c')
        stats_one = df_one['Cell_fraction'].describe()

        stats[key].append(c)
        stats['n'].append(df_one.shape[0])
        stats['min_freq'].append(stats_one['min'])
        stats['max_freq'].append(stats_one['max'])
        stats['mean_freq'].append(stats_one['mean'])
        stats['median_freq'].append(stats_one['50%'])
        stats['std_freq'].append(stats_one['std'])
        stats['SH'].append(SH(df_one['Cell_fraction']))
        stats['EI'].append(EI(df_one['Cell_fraction']))

    return pd.DataFrame(stats).set_index(key)


##


def SH(p):
    '''Shannon entropy of some pdist p.'''
    return - np.sum( np.log10(p) * p )


##


def EI(p):
    '''Expansion Index of some pdist p.'''
    return 1 -  ( SH(p) / np.log10(p.size) )


## 


def common(df):
    '''Median and mean common clones across samples...'''

    samples = df['sample'].cat.categories
    n_samples = samples.size

    counts = np.zeros((n_samples, n_samples))
    common_idxs = {}

    for i, x in enumerate(samples):
        GBC_x = set(df.query('sample == @x').index.to_list())
        for j, y in enumerate(samples):
            GBC_y = set(df.query('sample == @y').index.to_list())
            inter =  GBC_x & GBC_y
            common_idxs[f'{x}|{y}'] = inter
            counts[i, j] = len(inter)

    d = { 'median' : np.median(counts), 'mean' : np.mean(counts) }

    return d, counts, common_idxs
    

##


def common_in(df, n, normalize=False):
    '''% of clones shared by at least n samples.'''
    GBCs = df.index.unique()
    common_counts = []

    for x in GBCs:
        d = df.loc[x, 'sample']
        try:
            common_counts.append(df.loc[x, 'sample'].shape[0])
        except:
            common_counts.append(1)

    if normalize:
        return np.sum(np.array(common_counts) >= n) / GBCs.size
    else:
        return np.sum(np.array(common_counts) >= n) 


#


########################################################################

# Input data and paths
path_main = '/Users/IEO5505/Desktop/MDA_chemo_repro/clonal/'
path_data = path_main + '/data/'
path_results = path_main + '/results_and_plots/'

d = pd.read_excel(path_data + 'GBCs_freqs.xlsx', sheet_name=None)

# Format
for k in d:
    try:
        d[k] = d[k].loc[:, ['GBC', 'Cell_fraction']].assign(sample=k).dropna()
    except:
        print(d[k].columns)

df = pd.concat([ v for v in d.values() ], axis=0).set_index('GBC')
df = format_df(df)

# Colors:
cols = sns.color_palette('dark', n_colors=3)
colors = { k: v for k, v in zip(['in_vitro', 'treated', 'untreated'], cols)}

########################################################################

# Plots and numbers 
df['sample'].cat.categories.size # 21 samples
df['condition'].cat.categories.size # 2 conditions
df.index.unique().size # 15220 unique total clones from a single infection


##


# Fig n clones
#======================================================================#

fig, axs = plt.subplots(1, 2, figsize=(8, 5), sharey=True)

df_ = df.groupby('sample').size().to_frame(name='n')
bar(df_, 'n', c='grey', s=0.8, a=0.8, ax=axs[0], annot_size=5)
format_ax(df_, axs[0], title='Clones by sample', ylabel='n', 
    xticks=df_.index, rotx=90, log=True
)

df_ = df.groupby('condition').size().to_frame(name='n').reset_index()
bar(df_, 'n', by='condition', c=colors, s=0.8, a=0.8, ax=axs[1], annot_size=10)
format_ax(df_, axs[1], title='Clones by condition', 
    xticks=df_['condition'].unique(), rotx=90, log=True
)

fig.tight_layout()
fig.savefig(path_results + 'n_clones.pdf')

#======================================================================#


##


# Summaries
sample_df = stats_summary(df, 'sample')
sample_df.to_excel(path_results + 'summary_samples.xlsx')
condition_df = stats_summary(df, 'condition')
condition_df.to_excel(path_results + 'summary_conditions.xlsx')

d, counts, common_idxs = common(df) # 34.0, 76.93 median/mean
d, counts, common_idxs = common(df.query('condition == "untreated"')) # 0 26.47 median/mean
d, counts, common_idxs = common(df.query('condition == "treated"')) # 0 9.22 median/mean


##


# Fig n_commons
#======================================================================#

# All samples
counts = { 
    n : 
    common_in(df, n, normalize=False) \
    for n in range(2,22) # n samples
}

df_ = pd.DataFrame(data=counts.values(), index=counts.keys(), columns=['n'])

fig, ax = plt.subplots(figsize=(6, 5))
bar(df_, 'n', c='grey', s=0.8, a=0.8, ax=ax)
format_ax(df_, ax, title='n common clones', ylabel='n common clones', 
    xlabel='n tumors', xticks=df_.index
)
ax.set_ylim((0, 1900))
fig.savefig(path_results + 'n_commons.pdf')

#======================================================================#


##


# Fig difference in common clones distribution
#======================================================================#

# All samples
d, counts, common_idxs = common(df)

# Tests
combos = [ x.split('|') for x in common_idxs.keys() ]
counts_with_labels = { k : len(common_idxs[k]) for k in common_idxs }

df_ = pd.DataFrame(data=counts_with_labels.values(), columns=['n']).assign(
    between=counts_with_labels.keys()
)

test_tr = lambda x, y: bool(re.search('C1|C3', x)) & bool(re.search('C1|C3', y)) 
test_un = lambda x, y: bool(re.search('C2|C4', x)) & bool(re.search('C2|C4', y)) 

from itertools import starmap
test_tr_df = np.array(list(starmap(test_tr, [ x.split('|') for x in counts_with_labels.keys() ])))
test_un_df = np.array(list(starmap(test_un, [ x.split('|') for x in counts_with_labels.keys() ])))

df_ = df_.assign(condition=np.select([ test_tr_df, test_un_df ], ['treated', 'untreated']))
df_ = df_[df_['condition'] != '0']
c = { k : v for k, v in colors.items() if k != 'in_vitro' }

fig, ax = plt.subplots(figsize=(5, 5))
sns.stripplot(data=df_, x='condition', y='n', ax=ax, palette=c.values())
format_ax(df_, ax, title='Common clones by condition', 
    ylabel='n common clones', xlabel=''
)
add_wilcox(df_, 'condition', 'n', [('treated', 'untreated')], ax)
fig.savefig(path_results + 'n_commons_per_condition.pdf')

#======================================================================#



##


# Fig bubbles
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

df_ = df.reset_index().sample(df.shape[0])
df_['area_plot'] = rescale(df_['Cell_fraction']) * (3000-10) + 10

# Fig 
fig, ax = plt.subplots(figsize=(7, 7))

scatter(df_, 'GBC', 'sample', by='GBC', c=clones_colors, s='area_plot', a=0.5, ax=ax)
format_ax(df_, ax, title='Clones by sample', xlabel='Clones', xticks='')

fig.tight_layout()
fig.savefig(path_results + 'bubbles.pdf')

#======================================================================#


##


# Fig bubbles: 2. Only untreated PTs
#======================================================================#

import random

df_ = df.query('condition == "untreated"') 
n_clones = df_.index.unique().size

clones_colors = { 
    clone : color for clone, color in \
    zip(
        df_.index.unique(), 
        list(
            ''.join( ['#'] + [random.choice('ABCDEF0123456789') for i in range(6)] )  \
            for _ in range(n_clones)
        )
    )
}

df_ = df_.reset_index().sample(df_.shape[0])
df_['area_plot'] = rescale(df_['Cell_fraction']) * (3000-10) + 10

# Fig 
fig, ax = plt.subplots(figsize=(7, 7))

scatter(df_, 'GBC', 'sample', by='GBC', c=clones_colors, s='area_plot', a=0.5, ax=ax)
format_ax(df_, ax, title='Clones by sample', xlabel='Clones', xticks='')

fig.tight_layout()
fig.savefig('/Users/IEO5505/Desktop/revision_mini_NR_AC/results_and_plots/figures_final/bubbles_albi.pdf')

#======================================================================#


##


# Fig clustermap correlations
#======================================================================#

# Compute correlations
df_corr = df.reset_index().pivot_table(
    columns='sample', values='Cell_fraction', index='GBC'
)
df_corr.columns = list(df_corr.columns)
X = df_corr.values

C = np.zeros((df_corr.shape[1], df_corr.shape[1]))

for i, x in enumerate(df_corr.columns):
    for j, y in enumerate(df_corr.columns):
        clones_x = set(df_corr[x][~df_corr[x].isna()].index)
        clones_y = set(df_corr[y][~df_corr[y].isna()].index)
        intersection = list(clones_x & clones_y)
        freq_x = df_corr.loc[intersection, x]
        freq_y = df_corr.loc[intersection, y]
        C[i, j] = np.corrcoef(freq_x, freq_y)[0,1]

C.mean() # 0.24!!

df_ = pd.DataFrame(data=C, index=df_corr.columns, columns=df_corr.columns)

# Reorder
samples_order = [ 
        
    'Reference',
    'C1_I', 'C1_II', 'C1_III', 'C1_Null', 'C1_striped', 
    'C3_I', 'C3_II', 'C3_III', 'C3_striped', 
    'C2_I', 'C2_II', 'C2_III', 'C2_Null', 'C2_blu', 'C2_striped',
    'C4_I', 'C4_II', 'C4_III', 'C4_Null', 'C4_striped'
    
]

len(samples_order)
df_ = df_.reindex(samples_order)

# Tests
test_tr =  df_.index.str.startswith('C1')  | df_.index.str.startswith('C3')
test_un = df_.index.str.startswith('C2')  | df_.index.str.startswith('C4')

df_.loc[test_un, test_un].values.mean() # 0.46
df_.loc[test_tr, test_tr].values.mean() # 0.28

# Colors
f = [ df_.index == 'Reference', test_tr, test_un ]
v = ['in_vitro', 'treated', 'untreated']
row_colors = np.select(f, v).tolist()
row_colors = [ colors[k] for k in row_colors ]

# Round
df_ = np.round(df_, 2)

# Fig, 1
fig = plot_clustermap(df_, row_colors, palette='mako', 
    title='Clonal frequency correlation', label='Pearson correlation', 
    no_cluster=True, figsize=(8, 6.5), annot=True, annot_size=5, colors_ratio=0.02,
)
handles = create_handles(v, colors=colors.values())
fig.fig.subplots_adjust(left=0.2)
fig.fig.legend(handles, v, loc='lower center', bbox_to_anchor=(0.16, 0.5), ncol=1, frameon=False)
fig.ax_cbar.set_position((0.325, 0.05, 0.5, 0.02))
fig.savefig(path_results + 'corr_1.pdf')

# Fig, 2
fig = plot_clustermap(df_, row_colors, palette='mako', 
    title='Clonal frequency correlation', label='Pearson correlation', 
    no_cluster=False, figsize=(8.2, 6.5), annot=True, annot_size=5, colors_ratio=0.02,
)
handles = create_handles(v, colors=colors.values())
fig.fig.subplots_adjust(left=0.2)
fig.fig.legend(handles, v, loc='lower center', bbox_to_anchor=(0.12, 0.5), ncol=1, frameon=False)
fig.ax_cbar.set_position((0.325, 0.05, 0.5, 0.02))
fig.savefig(path_results + 'corr_2.pdf')

#======================================================================#


##


# Per sample metrics 
#======================================================================#

df_ = format_df(sample_df.reset_index())
df_ = df_.query('condition != "in_vitro"')
df_['condition'] = df_['condition'].cat.remove_unused_categories()
c = { k : v for k, v in colors.items() if k != 'in_vitro' }
pairs = [('untreated', 'treated')]

# Fig
fig, axs = plt.subplots(2, 2, figsize=(7, 6), sharex=True)

box(df_, 'condition', 'std_freq', c=c, ax=axs[0, 0], with_stats=True, pairs=pairs)
box(df_, 'condition', 'median_freq', c=c, ax=axs[0, 1], with_stats=True, pairs=pairs)
box(df_, 'condition', 'SH', c=c, ax=axs[1, 0], with_stats=True, pairs=pairs)
box(df_, 'condition', 'EI', c=c, ax=axs[1, 1], with_stats=True, pairs=pairs)

fig.tight_layout()
fig.savefig(path_results + 'metrics.pdf')

#======================================================================#


##

# Distribution by condition
#======================================================================#

df_ = df.reset_index()
df_['Cell_fraction'] = -np.log10(df_['Cell_fraction'])

fig, ax = plt.subplots(figsize=(5, 5))

hist(df_, 'Cell_fraction', n=round(df.shape[0]/25), by='condition', c=colors, a=0.5, ax=ax)
format_ax(df_, ax, title='Clones frequency distribution', xlabel='-log10(f)', ylabel='n')
ax.set(xlim=(1.5, 5.5))
handles = create_handles(colors.keys(), marker='o', colors=colors.values(), size=10, width=0.5)
ax.legend(handles, colors.keys(), frameon=False, loc='upper left', bbox_to_anchor=(0, 1))

fig.tight_layout()
fig.savefig(path_results + 'distributions.pdf')

#======================================================================#


##


########################################################################


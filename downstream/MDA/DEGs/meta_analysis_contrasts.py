"""
Script to create a consensus ranking of dominant pro-metastatic clones 
DEGs.
"""

import os
import re
import pickle
import numpy as np
import pandas as pd
import textalloc as ta
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


#


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_results = os.path.join(path_main, 'results', 'contrasts', 'consensus')

# contrasts_name = 'promet_clones_10_cells.pickle'

# Read pickle
# with open(os.path.join(path_results, contrasts_name), 'rb') as f:
#     results = pickle.load(f)
    
# Create consensus
# for pattern in ['AC_AC']:#:, 'NT_AC', 'AC_NT', 'NTA_NT']:
#     fig, df = consensus_volcano(
#         results, pattern=pattern, figsize=(7,7), return_df=True)
#     fig.savefig(os.path.join(path_results, f'{pattern}_volcano.png'))
#     (
#         df
#         .drop(columns=['log_evidence', 'to_annotate_positive', 'to_annotate_negative'])
#         .to_excel(os.path.join(path_results, f'{pattern}.xlsx'))
#     ) 
#     print(df[:5])


##

# Harmonize pathways
pathways = list(set(pd.read_csv(os.path.join(path_results, 'NT_AC_annotation.csv'), index_col=0)['Annotation'].unique()) | \
    set(pd.read_csv(os.path.join(path_results, 'AC_AC_annotation.csv'), index_col=0)['Annotation'].unique()) | \
    set(pd.read_csv(os.path.join(path_results, 'AC_AC_dominant_annotation.csv'), index_col=0)['Annotation'].unique()))
colors = { k:v for k,v in zip(pathways, sc.pl.palettes.default_20)}


##


# NT_AC
df = pd.read_excel(os.path.join(path_results, 'NT_AC.xlsx'), index_col=0)
df_anno = pd.read_csv(os.path.join(path_results, 'NT_AC_annotation.csv'), index_col=0)
df = df.join(df_anno)
df.loc['MANCR', 'Annotation'] = 'lncRNA'

t = 0.5
df = (
    df
    .assign(to_annotate_positive=lambda x: (x['effect_size']>=t) & (x['evidence']<=0.1))
    .assign(to_annotate_negative=lambda x: (x['effect_size']<=-t) & (x['evidence']<=0.1))
)
df['log_evidence'] = -np.log10(df['evidence']+0.00000001)
pathways = df['Annotation'].loc[lambda x: ~x.isna()].unique().tolist()

##

# Viz
fig, ax = plt.subplots(figsize=(6, 5))

df['Color_RGB'] = df['Annotation'].map(lambda x: colors[x] if x in colors else '#575555')

test = ~((df['to_annotate_positive']) | (df['to_annotate_negative']))
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='lightgrey', s=1)
test = (df['to_annotate_positive']) & (df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='r', s=20)
test = (df['to_annotate_negative']) & (df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='b', s=20)
test = (df['to_annotate_positive']) & (~df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], 
        c=df.loc[test]['Color_RGB'].values, s=20)
test = (df['to_annotate_negative']) & (~df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], 
        c=df.loc[test]['Color_RGB'], s=20)

ta.allocate_text(
    fig, ax, 
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['effect_size'],
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['log_evidence'],
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']].index,
    x_scatter=df['effect_size'], y_scatter=df['log_evidence'], 
    linecolor='black', textsize=8, 
    max_distance=0.5, linewidth=0.5, nbr_candidates=100
)

ax.set_xlim((-1.2,1.2))
add_legend(
    colors={ k : colors[k] for k in colors if k in df['Annotation'].values },
    ax=ax, label='Pathways', 
    ncols=1,
    artists_size=8,
    ticks_size=8,
    label_size=10,
    loc='center left', bbox_to_anchor=(1,.5)
)

ax.axis('off')
format_ax(ax, title='Consensus Adj (n clones=3)')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'NT_AC_annotated_volcano.png'), dpi=500)


##


# AC_AC
df = pd.read_excel(os.path.join(path_results, 'AC_AC.xlsx'), index_col=0)
df_anno = pd.read_csv(os.path.join(path_results, 'AC_AC_annotation.csv'), index_col=0)
df = df.join(df_anno)

t = 0.5
df = (
    df
    .assign(to_annotate_positive=lambda x: (x['effect_size']>=t) & (x['evidence']<=0.1))
    .assign(to_annotate_negative=lambda x: (x['effect_size']<=-t) & (x['evidence']<=0.1))
)
df['log_evidence'] = -np.log10(df['evidence']+0.00000001)
pathways = df['Annotation'].loc[lambda x: ~x.isna()].unique().tolist()

##

# Viz
fig, ax = plt.subplots(figsize=(6, 5))

df['Color_RGB'] = df['Annotation'].map(lambda x: colors[x] if x in colors else '#575555')

test = ~((df['to_annotate_positive']) | (df['to_annotate_negative']))
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='lightgrey', s=1)
test = (df['to_annotate_positive']) & (df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='r', s=20)
test = (df['to_annotate_negative']) & (df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='b', s=20)
test = (df['to_annotate_positive']) & (~df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], 
        c=df.loc[test]['Color_RGB'].values, s=20)
test = (df['to_annotate_negative']) & (~df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], 
        c=df.loc[test]['Color_RGB'], s=20)

ta.allocate_text(
    fig, ax, 
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['effect_size'],
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['log_evidence'],
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']].index,
    x_scatter=df['effect_size'], y_scatter=df['log_evidence'], 
    linecolor='black', textsize=8, 
    max_distance=0.5, linewidth=0.5, nbr_candidates=100
)

ax.set_xlim((-1.2,1.2))
add_legend(
    colors={ k : colors[k] for k in colors if k in df['Annotation'].values },
    ax=ax, label='Pathways', 
    ncols=1,
    artists_size=8,
    ticks_size=8,
    label_size=10,
    loc='center left', bbox_to_anchor=(1,.5)
)

ax.axis('off')
format_ax(ax, title='Consensus NeoAdj + Adj (n clones=3)')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'AC_AC_annotated_volcano.png'), dpi=500)


##


# AC_AC_dominant
df = pd.read_csv(os.path.join(path_results, 'AC_AC_dominant_consensus_df.csv'), index_col=0)
df_anno = pd.read_csv(os.path.join(path_results, 'AC_AC_dominant_annotation.csv'), index_col=0)
df = df.join(df_anno)

t = .8
df = (
    df
    .assign(to_annotate_positive=lambda x: (x['effect_size']>=t) & (x['evidence']<=0.1))
    .assign(to_annotate_negative=lambda x: (x['effect_size']<=-t) & (x['evidence']<=0.1))
)
df['log_evidence'] = -np.log10(df['evidence']+0.00000001)
pathways = df['Annotation'].loc[lambda x: ~x.isna()].unique().tolist()

##

# Viz
fig, ax = plt.subplots(figsize=(6, 5))

df['Color_RGB'] = df['Annotation'].map(lambda x: colors[x] if x in colors else '#575555')

test = ~((df['to_annotate_positive']) | (df['to_annotate_negative']))
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='lightgrey', s=1)
test = (df['to_annotate_positive']) & (df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='r', s=20)
test = (df['to_annotate_negative']) & (df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], c='b', s=20)
test = (df['to_annotate_positive']) & (~df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], 
        c=df.loc[test]['Color_RGB'].values, s=20)
test = (df['to_annotate_negative']) & (~df['Annotation'].isna())
ax.scatter(x=df.loc[test]['effect_size'], y=df.loc[test]['log_evidence'], 
        c=df.loc[test]['Color_RGB'], s=20)

ta.allocate_text(
    fig, ax, 
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['effect_size'],
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']]['log_evidence'],
    df.loc[lambda x: x['to_annotate_positive'] | x['to_annotate_negative']].index,
    x_scatter=df['effect_size'], y_scatter=df['log_evidence'], 
    linecolor='black', textsize=8, 
    max_distance=0.5, linewidth=0.5, nbr_candidates=100
)

#ax.set_xlim((-1.5,1.5))
add_legend(
    colors={ k : colors[k] for k in colors if k in df['Annotation'].values },
    ax=ax, label='Pathways', 
    ncols=1,
    artists_size=8,
    ticks_size=8,
    label_size=10,
    loc='center left', bbox_to_anchor=(1,.5)
)

ax.axis('off')
format_ax(ax, title='Consensus NeoAdj + Adj dominant clones (n=3)')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'AC_AC_dominant_annotated_volcano.png'), dpi=500)


##




































# Dominant resistant promet 
L = []

name = 'dom_AC_AC_PT_1.csv'
df_ = (
    pd.read_csv(os.path.join(path_main, 'data', name), index_col=0)
    .query('comparison == "g0_vs_g1"')
)
df_.loc[:, 'comparison'] = '_'.join([name, 'vs_rest'])
L.append(df_)

name = 'dom_AC_AC_PT_2.csv'
df_ = (
    pd.read_csv(os.path.join(path_main, 'data', name), index_col=0)
    .query('comparison == "g0_vs_g1"')
)
df_.loc[:, 'comparison'] = '_'.join([name, 'vs_rest'])
L.append(df_)

name = 'dom_AC_NT_PT_2.csv'
df_ = (
    pd.read_csv(os.path.join(path_main, 'data', name), index_col=0)
    .query('comparison == "g0_vs_g1"')
)
df_.loc[:, 'comparison'] = '_'.join([name, 'vs_rest'])
L.append(df_)

##

# All clones 
fig, df_ = consensus_volcano(L=L, pattern='all dominant', t=0.7, xlim=(-1.3,1.3),
            figsize=(7,7), return_df=True)
fig.savefig(os.path.join(path_results, 'all_dominant_volcano.png'))
df_.to_csv(os.path.join(path_results, 'all_dominant_consensus_df.csv'))


##


# Only ACAC
fig, df_ = consensus_volcano(L=L[:-1], pattern='ACAC', t=.75, xlim=(-3.5,3.5),
            figsize=(7,7), return_df=True)
fig.savefig(os.path.join(path_results, 'ACAC_dominant_volcano.png'))
df_.to_csv(os.path.join(path_results, 'ACAC_dominant_consensus_df.csv'))




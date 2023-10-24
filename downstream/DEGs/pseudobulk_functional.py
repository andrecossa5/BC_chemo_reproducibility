"""
Script to perform GSEA and visualization of pseudobulk DE.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
from itertools import chain, product
from functools import reduce
from plotting_utils import *
from BC_chemo_utils.tests import fastGSEA
from BC_chemo_utils.utils import read_model_results, read_GSEA
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_results = os.path.join(path_main, 'results', 'MDA', 'pseudobulk')

# Read DE results
df_fixed = (
    read_model_results(path_results, 'fixed')
    [['logFC', 'FDR', 'contrast', 'regressed_cc']]
    .assign(model='fixed')
)
df_mixed = (
    read_model_results(path_results, 'random')
    .rename(columns={'adj.P.Val':'FDR'})
    [['logFC', 'FDR', 'contrast', 'regressed_cc']]
    .assign(model='mixed')
)

# Merge and save
df = pd.concat([ df_fixed, df_mixed ])
df.to_csv(os.path.join(path_results, 'results.csv'))

# Read
df = pd.read_csv(os.path.join(path_results, 'results.csv'), index_col=0)


##


# GSEA on all DEGs lists
combos = product(
    df['contrast'].unique(), 
    df['model'].unique(), 
    df['regressed_cc'].unique()
)
for c, m, cc in combos:
    df_ = fastGSEA(df.query('contrast==@c and model==@m and regressed_cc==@cc')['logFC'])
    df_.to_csv(
        os.path.join(
            path_results, 'GSEA', 'GO_biological_process_2021', f'{c}_{m}_{cc}.csv'
        )
    )

# Read GSEA
df_gsea = read_GSEA(path_results)


##


# Make some sense out of it!

# Top represented genes across all contrasts: within top5, FDR and logFC
feat = 'logFC'
pd.Series(chain.from_iterable(
    df
    .groupby(['contrast', 'model', 'regressed_cc'])
    .apply(lambda x: x.sort_values(feat, ascending=False).index.to_list()[:5])
)).value_counts().to_frame('n').to_csv(
    os.path.join(path_results, f'top_represented_genes_{feat}.csv')
)

# Over-represented terms and genes in GSEA

# Terms
feat = 'p_adjusted'
pd.Series(chain.from_iterable(
    df_gsea
    .groupby(['contrast', 'model', 'cc'])
    .apply(lambda x: x.sort_values(feat, ascending=True).index.to_list()[:5])
)).value_counts().to_frame('n').to_csv(
    os.path.join(path_results, f'top_represented_terms_{feat}.csv')
)

# Leading genes
pd.Series(chain.from_iterable(
    df_gsea
    .groupby(['contrast', 'model', 'cc'])
    .apply(lambda x: x.head(5)['Lead_genes'].map(lambda x: x.split(';')))
)).value_counts().to_frame('n').to_csv(
    os.path.join(path_results, f'top_represented_leading_genes.csv')
)

# Differences for cc and mixed effects --> df with top 50 exclusive genes per model-cc combo
feat = 'FDR'
n = 50
perc = 75
grouped = (
    df
    .groupby(['contrast', 'model', 'regressed_cc'])
    .apply(lambda x: 
        x.query(f'logFC>{np.percentile(x["logFC"], perc)}')
        .sort_values(feat, ascending=True).index.to_list()[:n]
    )
)                                                                                                                                                                                                                                                                                                

d = {}
for contrast in grouped.index.levels[0].unique():
    for k in grouped.loc[contrast].index:
        diff = set(grouped.loc[contrast][k]) - \
        set(list(chain.from_iterable(grouped.loc[contrast].loc[lambda x: x.index!=k])))
        d[f'{contrast}_{k[0]}_{k[1]}'] = ';'.join(list(diff))

n_col = f'distinct_genes_in_top{n}(logFC>={perc} percentile), ranked by FDR'
df_distinct = pd.Series(d).to_frame(n_col)
df_distinct['perc_distinct'] = df_distinct[n_col].map(lambda x: len(x.split(';')) / n)
df_distinct = df_distinct.sort_values('perc_distinct', ascending=False)
df_distinct.to_csv(os.path.join(path_results, f'differences_in_top{n}.csv'))

# Differences for cc and mixed effects --> df with top 10 exclusive gsea terms per model-cc combo
feat = 'NES'
n = 10
grouped = (
    df_gsea
    .groupby(['contrast', 'model', 'cc'])
    .apply(lambda x:
        x.sort_values(feat, ascending=True).index.to_list()[:n]
    )
)                                                                                                                                                                                                                                                                                                

d = {}
for contrast in grouped.index.levels[0].unique():
    for k in grouped.loc[contrast].index:
        diff = set(grouped.loc[contrast][k]) - \
        set(list(chain.from_iterable(grouped.loc[contrast].loc[lambda x: x.index!=k])))
        d[f'{contrast}_{k[0]}_{k[1]}'] = ';'.join(list(diff))

n_col = f'distinct_terms_in_top{n}, ranked by NES'
df_distinct = pd.Series(d).to_frame(n_col)
df_distinct['perc_distinct'] = df_distinct[n_col].map(lambda x: len(x.split(';')) / n)
df_distinct = df_distinct.sort_values('perc_distinct', ascending=False)
df_distinct.to_csv(os.path.join(path_results, f'differences_in_top{n}_gsea.csv'))

# Consistency for cc and mixed effects --> df with top 50 shared genes per model-cc combo
feat = 'FDR'
n = 50
perc = 75
grouped = (
    df
    .groupby(['contrast', 'model', 'regressed_cc'])
    .apply(lambda x: 
        x.query(f'logFC>{np.percentile(x["logFC"], perc)}')
        .sort_values(feat, ascending=True).index.to_list()[:n]
    )
)     

d = {}
for contrast in grouped.index.levels[0].unique():
    for k in grouped.loc[contrast].index:
        L = grouped.loc[contrast]
        d[f'{contrast}_{k[0]}_{k[1]}'] = ';'.join(list(reduce(lambda x,y: set(x) & set(y), L)))
n_col = f'overlapping_genes_in_top{n}(logFC>={perc} percentile), ranked by FDR'
df_overlap = pd.Series(d).to_frame(n_col)
df_overlap['perc_overlap'] = df_overlap[n_col].map(lambda x: len(x.split(';')) / n)
new_idx = df_overlap.index.map(lambda x: '_'.join(x.split('_')[:-2]))
df_overlap = df_overlap.reset_index(drop=True).set_index(new_idx).drop_duplicates()
df_overlap = df_overlap.sort_values('perc_overlap', ascending=False)
df_overlap.to_csv(os.path.join(path_results, f'overlap_in_top{n}.csv'))


##
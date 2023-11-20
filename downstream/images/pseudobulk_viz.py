"""
UMAP of pseudobulk clones.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import scanpy as sc
from matplotlib.gridspec import GridSpec
from plotting_utils._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'pseudobulk')

# Read clustered
clustered = sc.read(os.path.join(path_data, 'clustered.h5ad'))
clustered.obs['GBC_sample'] = clustered.obs['GBC'].astype(str) + \
                              clustered.obs['sample'].astype(str)
GBC_size = clustered.obs.groupby('GBC_sample').size().loc[lambda x: x>=10]
clustered.obs['mock'] = clustered.obs['GBC'].astype(str) + clustered.obs['sample'].astype(str)

# Read top represented
name_files = [
    'top_represented_genes_logFC.csv',
    'top_represented_genes_FDR.csv',
    'top_represented_leading_genes.csv',
    'top_represented_terms_NES.csv',
]
names = ['genes_logFC', 'genes_FDR', 'leading_genes', 'terms']
d = {}
for name, name_file in zip(names, name_files):
    d[name] = pd.read_csv(os.path.join(path_results, name_file), index_col=0)


##


# Viz 
n = 20
fig = plt.figure(figsize=(14,4.5))
g = GridSpec(1,5, figure=fig, width_ratios=[1,1,1,1,1])

ax = fig.add_subplot(g[0])
stem_plot(d['genes_logFC'].head(n), 'n', ax=ax)
format_ax(ax, title='Genes ranked by logFC', xlabel='n')
ax.spines[['right', 'left', 'top']].set_visible(False)

ax = fig.add_subplot(g[1])
stem_plot(d['genes_FDR'].head(n), 'n', ax=ax)
format_ax(ax, title='Genes ranked by FDR', xlabel='n')
ax.spines[['right', 'left', 'top']].set_visible(False)

ax = fig.add_subplot(g[2])
stem_plot(d['leading_genes'].head(n), 'n', ax=ax)
format_ax(ax, title='Leading genes', xlabel='n')
ax.spines[['right', 'left', 'top']].set_visible(False)

ax = fig.add_subplot(g[4])
df_ = d['terms']
df_.index = df_.index.map(lambda x: ' '.join(x.split(' ')[:-1]))
df_.index = df_.index.map(lambda x: x[:40]+'...' if len(x)>40 else x)

df_.index = df_.index.map(lambda x: ' '.join(x.split(' ')[:6]+['...']))
stem_plot(df_.head(n), 'n', ax=ax)
format_ax(ax, title='GSEA terms (GO 2021)', xlabel='n')
ax.spines[['right', 'left', 'top']].set_visible(False)

fig.subplots_adjust(top=.9, bottom=.1, left=.1, right=.9, wspace=.6)
fig.savefig(os.path.join(path_results, 'pseudobulk_overrepresented_genes.png'), dpi=300)


##


# Read DE
name_files = [
    'results_treated_vs_untreated_PTs_nUMIs_mito_perc_G2.M.csv',
    'results_double_treated_vs_untreated_lungs_nUMIs_mito_perc_G2.M.csv',
    'results_PT_treated_vs_untreated_lungs_nUMIs_mito_perc.csv'
]
names = ['PT, treated vs untreated', 'lung, double- vs un-treated', 'lung, PT- vs un-treated']

d = {}
for name, name_file in zip(names, name_files):
    d[name] = pd.read_csv(os.path.join(path_results, 'fixed', name_file), index_col=0)


##


# Viz
for k in d:
    fig = volcano(d[k], effect_size='logFC', evidence='FDR', 
            n=10, annotate=True, xlim=(-12, 12))
    fig.tight_layout()
    fig.suptitle(k)
    fig.savefig(os.path.join(path_results, 'volcano', f'{k}.png'), dpi=300)


##

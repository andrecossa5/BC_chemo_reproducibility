"""
Definition and association of metastatic potential.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scanpy as sc
from plotting_utils._utils import *
from plotting_utils._colors import *
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
from BC_chemo_utils.clonal_utils import *
from BC_chemo_utils.tests import *
from Cellula.preprocessing._pp import scanpy_score
from matplotlib.gridspec import GridSpec
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal')

# Read cells meta
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))

# Calculate stats
df = calculate_clonal_SH_and_prevalence(adata)

# Retrieve PCs (and their genes info)
with open(os.path.join(
    path_data, 'PCs_expanded_clones', 'GBC_dom_promet_vs_nonmet_ALL_PCs'), 'rb') as f:
    d = pickle.load(f)

# Extract down genes
genes = d['gs']['g0_vs_g1']['PC13'].stats.index[-10:].tolist() + \
        d['gs']['g0_vs_g1']['PC11'].stats.index[-10:].tolist() + \
        d['gs']['g0_vs_g1']['PC5'].stats.index[-10:].tolist() + \
        d['gs']['g0_vs_g1']['PC17'].stats.index[-10:].tolist() + \
        d['gs']['g0_vs_g1']['PC19'].stats.index[-10:].tolist() 
genes = [ 
    g for g in genes if g not in \
    ['KMT2E-AS1', 'IRS1', 'OASL', 'PHLDA1', 'AL118516.1', 'HILPDA', 'CH25H']
]

# Check they are up in our expanded, treated clones
# adata.obs['sig'] = scanpy_score(adata, genes)
# clones = ["CGAGGGGATGGACTTCCG", "AACTCGACGCCTTATCAG", "GTCCACAGCACCTGCTCG"]
# cells = adata.obs.query('origin == "PT" and GBC in @clones').index
# not_cells = adata.obs.query('condition == "PT, treated" and GBC not in @clones').index
# adata.obs.loc[cells, 'sig'].median()
# adata.obs.loc[not_cells, 'sig'].median()

# Update df with expression for these genes
for gene in genes:
    agg = []
    for i in range(df.shape[0]):
        GBC = df.index[i]
        dataset = df.iloc[i,0]
        cells = adata.obs.query('GBC == @GBC and dataset == @dataset and origin=="PT"').index
        agg.append(np.mean(adata[cells, gene].X.A.flatten()))
    df[gene] = agg


##


# Linear model(s)

# Remove outliers
min_ = np.percentile(df['met_potential'], 5)
max_ = np.percentile(df['met_potential'], 95)
df = df.query('met_potential>=@min_ and met_potential<=@max_')
df['met_potential'] = np.log10(df['met_potential'])

# 1. Shannon entropies
covs = ['SH_PT', 'SH_lung', 'diff_SH']
sh_L = simple_lm(df, cov=covs, n=100)

# 2. Cell states prevalences
covs = df.columns[7:20]
states_L = simple_lm(df, cov=covs, n=100)

# 3. Top distinguishing PCs in expanded clones contrasts: top loading genes
covs = np.unique(genes)
pcs_genes_L = simple_lm(df, cov=covs, n=100)


##


# Viz
fig = plt.figure(figsize=(12,5.5)) 
gs = GridSpec(1,4,figure=fig, width_ratios=[1,0.5,1.5,2])

ax = fig.add_subplot(gs[0])
df_ = sh_L[2].loc[sh_L[2].index[::-1]]
for j in range(df_.shape[0]):
    c = 'r' if df_.iloc[j,:].median() >=0 else 'b'
    ax.plot(df_.iloc[j,:].median(), j, f'{c}o', markersize=5)
    ax.errorbar(df_.iloc[j,:].median(), j, xerr=df_.iloc[j,:].std(), c=c)

t = f'Clonal complexity \n MAPE: {np.mean(sh_L[0]):.2f} (+-{np.std(sh_L[0]):.2f}) \n R2: {np.mean(sh_L[1]):.2f} (+-{np.std(sh_L[1]):.2f})'
format_ax(
    ax, title=t, yticks=df_.index, xlabel='Linear model weight', 
    reduce_spines=True, title_size=11
)

ax = fig.add_subplot(gs[2])
df_ = states_L[2].loc[states_L[2].index[::-1]]

for j in range(df_.shape[0]):
    c = 'r' if df_.iloc[j,:8].median() >=0 else 'b'
    ax.plot(df_.iloc[j,:8].median(), j, f'{c}o', markersize=5)
    ax.errorbar(df_.iloc[j,:8].median(), j, xerr=df_.iloc[j,:8].std(), c=c)
t = f'Cell states frequency \n MAPE: {np.mean(states_L[0]):.2f} (+-{np.std(states_L[0]):.2f}) \n R2: {np.mean(states_L[1]):.2f} (+-{np.std(states_L[1]):.2f})'
format_ax(
    ax, title=t, yticks=df_.index, xlabel='Linear model weight', 
    reduce_spines=True, title_size=11
)

ax = fig.add_subplot(gs[3])
df_ = pcs_genes_L[2].loc[pcs_genes_L[2].index[::-1]]
for j in range(df_.shape[0]):
    c = 'r' if df_.iloc[j,:].median() >=0 else 'b'
    ax.plot(df_.iloc[j,:].median(), j, f'{c}o', markersize=5)
    ax.errorbar(df_.iloc[j,:].median(), j, xerr=df_.iloc[j,:].std(), c=c)

t = f'Top PCs genes \n MAPE: {np.mean(pcs_genes_L[0]):.2f} (+-{np.std(pcs_genes_L[0]):.2f}) \n R2: {np.mean(pcs_genes_L[1]):.2f} (+-{np.std(pcs_genes_L[1]):.2f})'
format_ax(
    ax, title=t, yticks=df_.index, xlabel='Linear model weight', reduce_spines=True,
    title_size=11, yticks_size=6
)

fig.subplots_adjust(left=.1, right=.9, top=.83, bottom=.1, wspace=.5)
fig.suptitle('Linear associaton with metastatic potential')
fig.savefig(os.path.join(path_results, 'linear_models_metpotential.png'), dpi=300)


##








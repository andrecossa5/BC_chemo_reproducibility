"""
Clone calling MDA dataset.
"""

import os
import sys
import pickle
import pandas as pd
import numpy as np
from itertools import chain
from plotting_utils._plotting_base import *
sys.path.append('/Users/IEO5505/Desktop/MI_TO/mito_preprocessing/bin/sc_gbc')
from helpers import *


##


# Args
path_data = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/MDA'
path_clones = os.path.join(path_data, 'clonal_info')


##


# Reformat meta: new metadata
meta = pd.read_csv(os.path.join(path_data, 'cells_meta_orig.csv'), index_col=0)

# Add dataset
meta['dataset'] = (
    meta['sample']
    .map(lambda x: '_'.join(x.split('_')[:2]) + '_' + x.split('_')[-1] )
)
meta['origin'] = np.where(meta['sample'].str.contains('mets'), 'lung', 'PT')
L = [
    (meta['sample'].str.startswith('AC_AC')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('AC_AC')) & (meta['origin'] == 'lung'),
    (meta['sample'].str.startswith('NT_NT')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('NT_NT')) & (meta['origin'] == 'lung'),
    (meta['sample'].str.startswith('NT_AC')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('NT_AC')) & (meta['origin'] == 'lung')
]
c = [
    'PT, treated', 'lung, double-treated', 'PT, untreated', 'lung, untreated',
    'PT, untreated', 'lung, single-treated'
]
meta['condition'] = np.select(L, c)
meta['condition'] = pd.Categorical(
    meta['condition'], 
    categories=[
        'PT, treated', 'PT, untreated', 'lung, double-treated', 
        'lung, single-treated', 'lung, untreated',
    ]
)


##


# Check single-samples
sample = 'AC_NT_mets_3'

# Params
correction_type = 'reference-free'
coverage_treshold=20
umi_treshold=5
p_treshold=.5
ratio_to_most_abundant_treshold=.5


##


# Read counts as pickle
with open(os.path.join(path_data, 'clonal_info', f'{sample}_counts.pickle'), 'rb') as p:
    COUNTS = pickle.load(p)

# Filter only QCed cells
counts = COUNTS[correction_type]
counts = counts.loc[
    counts['CBC']
    .isin(meta.query('sample==@sample')
    .index.map(lambda x: x.split('_')[0]))
].copy()


##


# Filter UMIs
counts = mark_UMIs(counts, coverage_treshold=coverage_treshold, nbins=50)
fig, ax = plt.subplots(figsize=(5,5))
viz_UMIs(counts, by='status', ax=ax, nbins=50)
fig.tight_layout()
plt.show()

# Get combos
df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')

# Filtering CBC-GBC
M, _ = filter_and_pivot(
    df_combos, 
    umi_treshold=umi_treshold, 
    p_treshold=p_treshold,  
    ratio_to_most_abundant_treshold=ratio_to_most_abundant_treshold
)

# GBC sets checks
sets = get_clones(M)
sets.head(20)
GBC_set = list(chain.from_iterable(sets['GBC_set'].map(lambda x: x.split(';')).to_list()))
redundancy = 1-np.unique(GBC_set).size/len(GBC_set)
occurrences = pd.Series(GBC_set).value_counts().sort_values(ascending=False)
occurrences.median()

# Get 1-GBC CBCs
unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
filtered_M = M.loc[unique_cells]
clones_df = get_clones(filtered_M)
cells_df = (
    filtered_M
    .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
    .to_frame('GBC_set')
)

# Final clones checks
print(f'# Final clones (i.e., distinct populations of uniquely barcoded cells only) checks \n')
print(f'- n starting CBC (STARSolo): {COUNTS[correction_type]["CBC"].unique().size}\n')
print(f'- n starting CBC (QC cells): {meta.query("sample==@sample").shape[0]}\n')
print(f'- n uniquely barcoded cells: {cells_df.shape[0]}\n')
print(f'- n clones: {clones_df.shape[0]}\n')
print(f'- n clones>=10 cells: {clones_df["n cells"].loc[lambda x:x>=10].size}\n')

# Viz p_poisson vs nUMIs
fig, ax = plt.subplots(figsize=(5,5))
scatter(df_combos, 'umi', 'p', by='max_ratio', marker='o', s=10, vmin=.2, vmax=.8, ax=ax, c='Spectral_r')
format_ax(
    ax, title='p Poisson vs nUMIs, all CBC-GBC combinations', 
    xlabel='nUMIs', ylabel='p', reduce_spines=True
)
ax.axhline(y=p_treshold, color='k', linestyle='--')
ax.text(.2, .9, f'Total CBC-GBC combo: {df_combos.shape[0]}', transform=ax.transAxes)
n_filtered = df_combos.query('status=="supported"').shape[0]
ax.text(.2, .86, 
    f'n CBC-GBC combo retained: {n_filtered} ({n_filtered/df_combos.shape[0]*100:.2f}%)',
    transform=ax.transAxes
)
fig.tight_layout()
plt.show()


##


# Save chosen treshold for last pipeline run.


##

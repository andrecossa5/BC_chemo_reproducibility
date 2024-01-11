"""
Final cell_assignment MDA dataset.
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
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'MDA')
path_clones = os.path.join(path_data, 'clonal_info')
path_params = os.path.join(path_data, 'clone_assignment_params_small.csv')


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
    (meta['sample'].str.startswith('NT_AC')) & (meta['origin'] == 'lung'),
    (meta['sample'].str.startswith('AC_NT')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('AC_NT')) & (meta['origin'] == 'lung')
]
c = [
    'PT, treated', 'lung, double-treated', 'PT, untreated', 'lung, untreated',
    'PT, untreated', 'lung, single-treated', 'PT, treated', 'lung, single-treated'
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


# Get saved coverage_tresholds
df_params = pd.read_csv(path_params, index_col=0)

# Other params
correction_type = 'reference-free'
umi_treshold = 5
p_treshold = 1 
max_ratio_treshold = .75,
normalized_abundance_treshold = .75


# Get all cells, and visualization for all samples
CELL_DFs = []
dist_reads = plt.figure(figsize=(20,13))
dist_combos = plt.figure(figsize=(20,13))

# Read saved sample specific treshold
for i,sample in enumerate(df_params.index):

    # Get params
    d = df_params.loc[sample].to_dict()

    # Get counts
    with open(os.path.join(path_data, 'clonal_info', f'{sample}_counts.pickle'), 'rb') as p:
        COUNTS = pickle.load(p)

    # Filter only QCed cells
    counts = COUNTS['reference-free']
    counts = counts.loc[
        counts['CBC']
        .isin(meta.query('sample==@sample')
        .index.map(lambda x: x.split('_')[0]))
    ].copy()

    # Filter UMIs
    counts = mark_UMIs(counts, coverage_treshold=int(d['coverage_treshold']), nbins=50)

    # Add viz dist_reads
    ax = dist_reads.add_subplot(4,6,i+1)
    viz_UMIs(counts, by='status', ax=ax, nbins=50)
    ax.set(title=sample)

    # Get combos
    df_combos = get_combos(counts, gbc_col=f'GBC_{correction_type}')

    # Filtering CBC-GBC
    M, _ = filter_and_pivot(
        df_combos, 
        umi_treshold=umi_treshold, 
        p_treshold=p_treshold,  
        max_ratio_treshold=max_ratio_treshold,
        normalized_abundance_treshold=normalized_abundance_treshold
    )

    # Add viz dist_combos
    ax = dist_combos.add_subplot(4,6,i+1)
    sns.kdeplot(data=df_combos, x='normalized_abundance', y='max_ratio', ax=ax)
    ax.axvline(x=normalized_abundance_treshold, color='k')
    ax.axhline(y=max_ratio_treshold, color='k')
    ax.set(title=sample)

    # GBC sets checks
    sets = get_clones(M)
    unique_cells = (M>0).sum(axis=1).loc[lambda x: x==1].index
    ax.text(.6,.1, f'n cells: {len(unique_cells)}', transform=ax.transAxes)
    filtered_M = M.loc[unique_cells]
    clones_df = get_clones(filtered_M)
    cells_df = (
        filtered_M
        .apply(lambda x: filtered_M.columns[x>0][0], axis=1)
        .to_frame('GBC_set')

    )
    cells_df.index = cells_df.index.map(lambda x: x + f'_{sample}')

    # Append
    CELL_DFs.append(cells_df)


##
    

# Save figs and update meta
dist_reads.tight_layout()
dist_reads.savefig(os.path.join(path_main, 'results', 'MDA', 'clonal', 'UMI_tresholds.png'), dpi=200)
dist_combos.tight_layout()
dist_combos.savefig(os.path.join(path_main, 'results', 'MDA', 'clonal', 'filtered_combos.png'), dpi=200)

# New meta
meta_new = pd.concat(CELL_DFs).join(meta).rename(columns={'GBC_set':'GBC'})
meta_new['sample'].unique().size

# Add batches and reformat columns
batches_df = pd.read_csv(os.path.join(path_data, 'batches.csv'))
batches_df.columns = ['sample', 'seq_run', 'infection', 'exp_condition']
meta_new = meta_new.merge(batches_df[['sample', 'seq_run', 'infection']], on='sample')

# Inspect final metadata

# n cells by seq_run and infection
(
    meta_new
    .query('origin=="lung"')
    .groupby('seq_run')
    ['sample'].value_counts()
    .to_frame('n').reset_index()
    .groupby('seq_run')['n'].median()
)

# n clones by condition
(
    meta_new.groupby('sample')
    ['GBC'].nunique()
    .to_frame('nGBC')
    .reset_index()
    .merge(
        meta_new[['sample', 'condition']].drop_duplicates(), 
        how='left'
    )
    .groupby('condition')
    ['nGBC'].median()
)


##


# Save final metadata
meta_new.to_csv(os.path.join(path_data, 'cells_meta.csv'))


##
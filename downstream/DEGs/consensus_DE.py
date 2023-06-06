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
from BC_chemo_utils.plotting import *


#


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_results = os.path.join(path_main, 'results', 'contrasts')
contrasts_name = 'promet_clones_10_cells.pickle'

# Read pickle
with open(os.path.join(path_results, contrasts_name), 'rb') as f:
    results = pickle.load(f)
    
# Create consensus
for pattern in ['AC_AC']:#:, 'NT_AC', 'AC_NT', 'NTA_NT']:
    fig, df = consensus_volcano(
        results, pattern=pattern, figsize=(7,7), return_df=True)
    fig.savefig(os.path.join(path_results, f'{pattern}_volcano.png'))
    (
        df
        .drop(columns=['log_evidence', 'to_annotate_positive', 'to_annotate_negative'])
        .to_excel(os.path.join(path_results, f'{pattern}.xlsx'))
    ) 
    print(df[:5])


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




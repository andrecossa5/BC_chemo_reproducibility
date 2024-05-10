"""
Bulk clonal dynamics under treatment.
"""

import os
import pandas as pd
import numpy as np
from plotting_utils._plotting_base import *


##


# Args
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'CTCs')
path_results = os.path.join(path_main, 'results', 'CTCs', 'clonal')


##


# Read data
df = pd.read_csv(os.path.join(path_data, 'bulk_GBC_NT_AC.csv'), index_col=0)
df = df.loc[~(df['sample'].str.contains('Late') | df['sample'].str.contains('Early'))]

# Format data
tests = [ df['sample'].str.startswith('PT'), df['sample'].str.startswith('CTC'), df['sample'].str.startswith('LN') ]
options = ['PT', 'CTC', 'LN']
df['origin'] = pd.Categorical(np.select(tests, options), categories=options)
df['condition'] = np.where(df['sample'].str.contains('NT'), 'NT', 'AC')
df['dataset'] = df['sample'].map(lambda x: x.split('_')[1])
df['freq'] = df.groupby('sample')['read_count'].transform(lambda x: x/x.sum())
df = df.reset_index().rename(columns={'index':'GBC'})

# To wide
df_grouped = (
    df.groupby('dataset')
    .apply(
        lambda x: x.pivot_table(
            values='freq', columns='origin', index='GBC', fill_value=0, dropna=False)
    )
    .reset_index()
    .set_index('GBC')
    .assign(
        is_longitudinal=lambda x: (x.loc[:,options]>0).sum(axis=1)>1,
        CTC_PT=lambda x: x['CTC'] / (x['PT']+0.000001),
        CTC_LN=lambda x: x['CTC'] / (x['LN']+0.000001)
    )
    .reset_index()
    .merge(df[['GBC', 'dataset', 'condition']])
)

# Questions
df.groupby('sample')['GBC'].nunique().sort_values(ascending=False)
df_grouped.loc[df_grouped['PT']>0]['condition'].value_counts()                              # PT
df_grouped.loc[df_grouped['CTC']>0]['condition'].value_counts()                             # CTC
df_grouped.loc[df_grouped['LN']>0]['condition'].value_counts()                              # LN
df_grouped['condition'].value_counts()                                                      # longitudinal
df_grouped.loc[df_grouped['is_longitudinal']]['condition'].value_counts()                   # longitudinal
df_grouped.loc[np.sum(df_grouped[options]>0, axis=1)==3]['condition'].value_counts()        # all 3 sites

# AAAAAAAAADICOPORUEO
df_1 = df_grouped.loc[df_grouped['is_longitudinal']].groupby('dataset')['GBC'].nunique().sort_values(ascending=False)
df_2 = df_grouped.groupby('dataset')['GBC'].nunique().sort_values(ascending=False).loc[df_1.index]


# Save
df_grouped.to_csv(os.path.join(path_results, 'df_longitudinal_CTCs.csv'))


##
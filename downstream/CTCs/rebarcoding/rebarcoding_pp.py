"""
Rebarcoding
"""

import os
import numpy as np
import pandas as pd


##


# Path main
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data/CTCs/rebarcoding')
path_results = os.path.join(path_main, 'results/CTCs/rebarcoding')

# Read
first_barcoding = pd.read_csv(os.path.join(path_data, 'bulk_GBC_first_barcoding.csv'), index_col=0)
second_barcoding = pd.read_csv(os.path.join(path_data, 'bulk_GBC_rebarcoding.csv'), index_col=0)
ctc_original = first_barcoding.loc[first_barcoding['sample']=='CTC_NT3']
ctc_line = first_barcoding.loc[first_barcoding['sample']=='CTC_NT3_line']
ref_rebarcoding = second_barcoding.loc[second_barcoding['sample']=='re_barc_ref']
lungs_rebarcoding = second_barcoding.loc[second_barcoding['sample']!='re_barc_ref']
spikeins = pd.read_csv(os.path.join(path_data, 'spikein_16.12.24.csv'), header=None)[0].values

# Remove spike-ins re-calculate frequency
df = pd.concat([ctc_original, ctc_line, ref_rebarcoding, lungs_rebarcoding])
df = df.loc[~df.index.isin(spikeins)]
df['f'] = df.groupby('sample')['read_count'].transform(lambda x: x / x.sum())


##



df['sample'].unique()


GBCs_pre = set(df.query('sample=="CTC_NT3"').index) | set(df.query('sample=="CTC_NT3_line"').index)
GBC_re_barc_ref = set(df.query('sample=="re_barc_ref"').index)


len(GBCs_pre)
len(GBCs_pre - GBC_re_barc_ref)

df.query('sample=="CTC_NT3_line"').head()
df.query('sample=="re_barc_ref"').head()


df.query('sample=="CTC_NT3"').head(5)
df.query('sample=="CTC_NT3_line"').head(5)


df.query('sample=="re_barc_ref"').loc[df.query('sample=="CTC_NT3_line"').index[:5]]

df.query('sample=="CTC_NT3_line"').shape

df.query('sample=="re_barc_ref"').loc[df.query('sample=="re_barc_ref"').index.isin(df.query('sample=="CTC_NT3_line"').index)].shape

df.loc[df['sample'].str.contains('re_barc_lung')].groupby('sample').apply(lambda x: x.head(5))

df_ =df.loc[df['sample'].str.contains('re_barc_lung')].loc[lambda x: ~x.index.isin(GBCs_pre)]


df_.groupby('sample').apply()


from plotting_utils._plotting_base import *
bb_plot
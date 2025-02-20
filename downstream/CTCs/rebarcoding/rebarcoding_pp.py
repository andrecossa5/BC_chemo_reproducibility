"""
Rebarcoding
"""

import os
import re
import numpy as np
import pandas as pd
import random
import pickle
from plotting_utils._plotting import *
matplotlib.use('macOSX')

##


# Path main
path_main = '/Users/ieo7295/Desktop/BC_chemo_reproducibility'
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

# Remove barcode CTC_NT3 & CTC_NT3_line from re_barc and re-calculate frequency
df_nt=df[df['sample'].isin(['CTC_NT3','CTC_NT3_line'])].reset_index().rename(columns={'index':'GBC'})
df=df[~ df.index.isin(df_nt['GBC'])]
df['f'] = df.groupby('sample')['read_count'].transform(lambda x: x / x.sum())

df_freq=(df.reset_index().rename(columns={'index':'GBC'})  
        .groupby('sample')
        .apply(lambda x: x.assign(    
        cum_freq=(x['read_count'] / x['read_count'].sum()).cumsum()
    ))
    .reset_index(drop=True)
)


#bubble plot 
categories=['re_barc_lung_1', 're_barc_lung_2',
       're_barc_lung_3', 're_barc_lung_4', 're_barc_lung_5',
       're_barc_lung_6', 're_barc_lung_7', 're_barc_lung_8',
       're_barc_lung_9', 're_barc_lung_10', 're_barc_lung_11','re_barc_ref']

df_freq['sample'] = pd.Categorical(df_freq['sample'], categories=categories)
df_freq.sort_values(by=['sample'], inplace=True)

#Random colors for clones
# clones = df_freq['GBC'].unique()
# random.seed(1235)
# clones_colors = { 
#     clone : color for clone, color in \
#     zip(
#         clones, 
#         list(
#             ''.join( ['#'] + [random.choice('ABCDEF0123456789') for i in range(6)] )  \
#             for _ in range(clones.size)
#         )
#     )
# }
# with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'wb') as f:
#     pickle.dump(clones_colors, f)

with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'rb') as f:
    clones_colors = pickle.load(f)

df_freq['area_plot'] = df_freq['f'] * (3000-5) + 5

fig, ax = plt.subplots(figsize=(6, 6))
scatter(df_freq, 'GBC', 'sample', by='GBC', c=clones_colors, s='area_plot', a=0.5, ax=ax)
format_ax(ax, title='Clones by sample', xlabel='Clones', xticks='')
ax.text(.3, .23, f'n clones total: {df_freq["GBC"].unique().size}', transform=ax.transAxes)
fig.tight_layout()
fig.savefig(os.path.join(path_results,'bubble_plot.png'),dpi=300)


#sample, Shannon entropy, origin, n_clones
SH = []
for s in df_freq['sample'].unique():
    df_ = df_freq.query('sample==@s')
    x = df_['f']
    SH.append(-np.sum( np.log10(x) * x ))

df_sample = (
    pd.Series(SH, index=df_freq['sample'].unique())
    .to_frame('SH')
    .sort_values(by='SH', ascending=False)
    .reset_index().rename(columns={'index':'sample'})
    .merge(df_freq[['sample']], on='sample')
    .drop_duplicates()
    .set_index('sample')
    .assign(
        n_clones=lambda df_: df_.index.map(
            lambda s: df_freq[df_freq['sample'] == s].index.nunique()
        ))
)


#bar plot n_clones by sample 
df_sample_sorted= df_sample.sort_values(by='n_clones', ascending=True)
fig, ax = plt.subplots(figsize=(10,4.5))
bar(df_sample_sorted, 'n_clones', 'sample', s=.70, c='k', a=.7, ax=ax)
format_ax(ax=ax, title='n clones by sample', ylabel='n_clones', xticks=df_sample_sorted.index, rotx=90)
ax.spines[['left', 'top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'n_clones.png'), dpi=500)

#bar_plot of SH
fig, ax = plt.subplots(figsize=(10,4.5))
bar(df_sample, 'SH', 'sample', s=.70, c='k', a=.7, ax=ax)
format_ax(ax=ax, title='SH by sample', ylabel='SH', xticks=df_sample.index, rotx=90)
ax.spines[['left', 'top', 'right']].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'barplot_sh.png'), dpi=500)


#box,strip n_clones by condition
df_sample=df_sample.reset_index()
df_sample['condition'] = df_sample['sample'].apply(lambda x: 're_barc_lung' if 're_barc_lung_' in x else 're_barc_ref')
fig, ax = plt.subplots(figsize=(6,6))
box(df_sample, x='condition', y='n_clones', ax=ax, with_stats=False, 
    #order=['re_barc_ref','re_barc_lung_1', 're_barc_lung_2',
       #'re_barc_lung_3', 're_barc_lung_4', 're_barc_lung_5',
       #'re_barc_lung_6', 're_barc_lung_7', 're_barc_lung_8',
       #'re_barc_lung_9', 're_barc_lung_10', 're_barc_lung_11']
    order=['re_barc_ref','re_barc_lung']
)
strip(df_sample, x='condition', y='n_clones', ax=ax, order=['re_barc_ref','re_barc_lung']
    #   order=['re_barc_ref','re_barc_lung_1', 're_barc_lung_2',
    #    're_barc_lung_3', 're_barc_lung_4', 're_barc_lung_5',
    #    're_barc_lung_6', 're_barc_lung_7', 're_barc_lung_8',
    #    're_barc_lung_9', 're_barc_lung_10', 're_barc_lung_11']
    , c='k')

format_ax(ax=ax, title='n_clones by condition', ylabel='n_clones', rotx=90, reduce_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'n_clones_condition_alt.png'), dpi=300)

#box,strip SH by condition
fig, ax = plt.subplots(figsize=(6,6))
box(df_sample, x='condition', y='SH', ax=ax, with_stats=False, order=['re_barc_ref','re_barc_lung'] 
    #order=['re_barc_ref','re_barc_lung_1', 're_barc_lung_2',
    #    're_barc_lung_3', 're_barc_lung_4', 're_barc_lung_5',
    #    're_barc_lung_6', 're_barc_lung_7', 're_barc_lung_8',
    #    're_barc_lung_9', 're_barc_lung_10', 're_barc_lung_11']
)
strip(df_sample, x='condition', y='SH', ax=ax, order=['re_barc_ref','re_barc_lung']
#order=['re_barc_ref','re_barc_lung_1', 're_barc_lung_2',
#        're_barc_lung_3', 're_barc_lung_4', 're_barc_lung_5',
#        're_barc_lung_6', 're_barc_lung_7', 're_barc_lung_8',
#        're_barc_lung_9', 're_barc_lung_10', 're_barc_lung_11']
     , c='k')
format_ax(ax=ax, title='Shannon Entropy samples', ylabel='SH', rotx=90, reduce_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'SH_alt.png'), dpi=300)



#heatmap of common clones 
df=df.reset_index().rename(columns={'index':'GBC'})
d = {sample: set(df[df['sample'] == sample]['GBC']) for sample in df['sample'].unique()}

n=len(d)
C=np.zeros((n,n))

sample=list(d.keys())
for i,x in enumerate(sample):
    for j,y in enumerate(sample):
        if i >= j:
            common_clones = len(d[x] & d[y])
            C[i,j] = C[j,i] = common_clones

df_cc=pd.DataFrame(C, index=sample, columns=sample) 

order = [
    're_barc_ref', 're_barc_lung_1', 're_barc_lung_2', 're_barc_lung_3',
    're_barc_lung_4', 're_barc_lung_5', 're_barc_lung_6', 're_barc_lung_7',
    're_barc_lung_8', 're_barc_lung_9', 're_barc_lung_10', 're_barc_lung_11'
]

prefix_pattern = re.compile(r'^(%s)' % '|'.join(order))
ordered_samples = sorted(
    {s for s in df_cc.index if prefix_pattern.match(s)},
    key=lambda x: (
        order.index(prefix_pattern.match(x).group(1)),
        int(re.search(r'_(\d+)$', x).group(1)) if re.search(r'_(\d+)$', x) else float('inf') 
    )
)
df_reordered = df_cc.loc[ordered_samples, ordered_samples]
vmin, vmax= 0,2500
fig, ax = plt.subplots(figsize=(10,8))
plot_heatmap(df_reordered, ax=ax, annot=True, title='n common clones', x_names_size=8, y_names_size=8, annot_size=5.5, cb=False)
sns.heatmap(data=df_reordered, ax=ax, robust=True, cmap="mako", vmin=vmin, vmax=vmax, fmt='.2f',cbar_kws={'fraction':0.05, 'aspect':35, 'pad': 0.02})
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'common_rt.png'), dpi=300)
plt.show()




























GBC_lung_11= set(df.query('sample=="re_barc_lung_11"').index) 
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
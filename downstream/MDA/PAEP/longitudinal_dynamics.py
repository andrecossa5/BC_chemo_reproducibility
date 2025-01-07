"""
PAEP into longitudinal clones
"""

import os
import pandas as pd
import numpy as np
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Read
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA')


##


# Read
paep = pd.read_csv(os.path.join(path_data, 'paep.csv'), index_col=0)
meta = pd.read_csv(os.path.join(path_data, 'full_embs.csv'), index_col=0)
longitudinal = pd.read_csv(os.path.join(path_data, 'longitudinal_clones.csv'), index_col=0)
longitudinal = pd.merge(
    longitudinal.melt(
        id_vars=['dataset', 'GBC'], value_vars=['PT_count', 'lung_count'], 
        var_name='origin', value_name='n cell'
    ).assign(origin=lambda x: x['origin'].map(lambda x: x.split('_')[0])),
    longitudinal.melt(
        id_vars=['dataset', 'GBC'], value_vars=['PT_freq', 'lung_freq'], 
        var_name='origin', value_name='freq'
    ).assign(origin=lambda x: x['origin'].map(lambda x: x.split('_')[0]))
    , on=['dataset', 'GBC', 'origin']
)
meta['branch'] = meta['dataset'].map(lambda x: '_'.join(x.split('_')[:-1]))
meta['PAEP'] = paep['expr']


##


# Viz

# Per branch
fig, ax = plt.subplots(figsize=(5,5))
sns.boxplot(data=meta.query('branch!="AC_NT"'), x='branch', y='PAEP', hue='origin', ax=ax)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'PAEP_per_sample_and_branch.png'), dpi=500)


##


# Per clone
df_ = longitudinal.merge(
    meta.groupby(['GBC', 'dataset', 'origin'])['PAEP'].median().reset_index(), 
    on=['GBC', 'dataset', 'origin']
)
df_['branch'] = df_['dataset'].map(lambda x: '_'.join(x.split('_')[:-1]))
df_['origin'] = np.where(df_['origin']=='PT',1,2)

# Per branch, across longitudinal clones
fig, axs = plt.subplots(1,3,figsize=(5,4), sharey=True)

for gbc, group in df_.query('branch=="NT_NT"').groupby("GBC"):
    ax = axs[0]
    ax.plot(group["origin"], group["PAEP"], marker='o', linestyle='-', color='k')
    ax.set_xticks([1, 2], labels=['PT', 'lung'])
    ax.set_xlim((0.5, 2.5))
    ax.spines[['right', 'top']].set_visible(False)
    ax.set(title='PT- and lung- \nuntreated', ylabel='PAEP expression')
for gbc, group in df_.query('branch=="NT_AC"').groupby("GBC"):
    ax = axs[1]
    ax.plot(group["origin"], group["PAEP"], marker='o', linestyle='-', color='k')
    ax.set_xticks([1, 2], labels=['PT', 'lung'])
    ax.set_xlim((0.5, 2.5))
    ax.spines[['right', 'top']].set_visible(False)
    ax.set(title='PT- untreated and \nlung- treated')
for gbc, group in df_.query('branch=="AC_AC"').groupby("GBC"):
    ax = axs[2]
    ax.plot(group["origin"], group["PAEP"], marker='o', linestyle='-', color='k')
    ax.set_xticks([1, 2], labels=['PT', 'lung'])
    ax.set_xlim((0.5, 2.5))
    ax.spines[['right', 'top']].set_visible(False)
    ax.set(title='PT- and met- \ntreated')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'PAEP_across_logitudinal_clones.png'), dpi=500)


##


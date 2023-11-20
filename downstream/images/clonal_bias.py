"""
Clonal bias.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.gridspec import GridSpec
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
from BC_chemo_utils.clonal_utils import SH
from BC_chemo_utils.tests import compute_enrichment
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal')

# Read adata and format degs for dotplot
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
adata.obs = adata.obs.rename(columns={'cell_states':'cell_state'})
meta = adata.obs

# Take out only clone-sample combos with > 10 cells
df_good_clones = (
    meta
    .groupby(['dataset', 'origin', 'GBC'])
    .size().loc[lambda x:x>=10]
    .reset_index(name='ncells')
)

L = []
for i in range(df_good_clones.shape[0]):

    dataset, origin, GBC = df_good_clones.iloc[i,:3].values
    s = np.where( 
        (meta['dataset']==dataset) & (meta['origin']==origin) & (meta['GBC']==GBC),
        'g', 'other'
    )

    df_ = pd.DataFrame({'group' : pd.Categorical(s), 'cell_state' : meta['cell_state']})
    L.append(
        compute_enrichment(df_, 'cell_state', 'group', 'g')
        [['odds_ratio', 'group']].assign(clonal_pop=f'{dataset}_{origin}_{GBC}')
    )

df_ = pd.concat(L)
df_ = df_.query('odds_ratio<=50')
df_ = (
    df_.groupby('group')['odds_ratio']
    .describe()[['50%', 'std']]
    .sort_values('50%')
)


##


# SH distribution
good_clones = meta.groupby('GBC').size().loc[lambda x:x>=10].index
meta_ = meta.query('GBC in @good_clones')
cross = pd.crosstab(meta_['GBC'], meta_['cell_state'], normalize=0)
sh = pd.Series(
    [ SH(cross.iloc[i,:].values+.0000001) for i in range(cross.shape[0]) ],
    index=cross.index
).to_frame('SH')

max_theorical = SH([1/cross.shape[1] for _ in range(cross.shape[1])])


##


fig = plt.figure(figsize=(12,5)) 
gs = GridSpec(1,2,figure=fig,width_ratios=[1,2])

ax = fig.add_subplot(gs[0])
for j in range(df_.shape[0]):
    ax.plot(df_.iloc[j,0], j, 'ko', markersize=10)
    ax.errorbar(df_.iloc[j,0], j, xerr=df_.iloc[j,1], c='k')

format_ax(
    ax, title='Cell state clonal bias', yticks=df_.index, 
    ylabel='Cell states', xlabel='odds ratio', reduce_spines=True
)

ax = fig.add_subplot(gs[1])
sns.kdeplot(data=sh, x='SH', ax=ax, fill=True)
format_ax(
    ax, title='Shannon Entropy of cell_state distributions', xlabel='SH', 
    ylabel='Density', reduce_spines=True
)
ax.vlines(x=max_theorical, ymin=0, ymax=4, linestyles='dashed')
fig.tight_layout()

fig.savefig(os.path.join(path_results, 'clonal_bias.png'), dpi=500)


##
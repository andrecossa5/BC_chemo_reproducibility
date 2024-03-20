"""
Longitudinal dynamics CTCs.
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import norm
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


def gaussian_smooth(x, y, grid, sd):
    weights = np.transpose([norm.pdf(grid, m, sd) for m in x])
    weights = weights / weights.sum(0)
    return (weights * y).sum(1)


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'CTCs')
path_results = os.path.join(path_main, 'results', 'CTCs')

# Read data
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
df_ = (
    meta.groupby('sample')
    .apply(lambda x: x['GBC'].value_counts(normalize=True))
    .reset_index()
    .pivot_table(columns='sample', index='level_1', values='GBC', fill_value=0)
)
# meta[['dataset', 'origin']].value_counts()

# Fishplot, 3 timepoints long sample
d = {
    'Dataset 1, ealy' : ['PT_1_early', 'CTC_1_early'],
    'Dataset 1, late' : ['PT_1_late', 'CTC_1_late', 'lung_1_late'],
    'Dataset 2, late' : ['PT_2_late', 'CTC_2_late', 'lung_2_late']
}

# Viz
for dataset, sample_list in d.items():

    df_long = df_[sample_list]
    df_long.sort_values(sample_list[-1], ascending=False)

    colors = { k:v for k,v in zip(df_long.index, sc.pl.palettes.godsnot_102[::-1])}
    x = np.arange(df_long.shape[1])
    y = [ np.array(x) for x in df_long.values.tolist() ]
    grid = np.linspace(-1.3, 4, num=500)
    y_smoothed = [ gaussian_smooth(x, y_, grid, .35) for y_ in y ]

    fig, ax = plt.subplots(figsize=(10,5))
    ax.stackplot(grid, y_smoothed, baseline="sym", colors=colors.values())
    format_ax(ax, xticks=sample_list, ylabel='Prevalence')
    n_total = np.sum((df_long>0).sum(axis=1)>0)
    n_long = np.sum((df_long>0).sum(axis=1)==df_long.shape[1])
    ax.set_title(
        f'{dataset}: n clones {n_total}, n longitudinal: {n_long}', loc='left')
    ax.spines[['right', 'bottom', 'top']].set_visible(False)
    for l in x:    
        ax.axvline(x=l, color='k', linewidth=.5, linestyle='dashed')
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, f'{dataset}_dynamics.png'), dpi=400)


##
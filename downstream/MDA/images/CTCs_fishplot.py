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
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA','CTCs')

# Read data
df = pd.read_csv(os.path.join(path_data, 'CTCs', 'CTCs_clonal_prevalences.csv'), index_col=0)
df_ = (
    df[['cellular_prevalence_wi', 'sample']]
    .pivot(columns='sample', values='cellular_prevalence_wi')  
) 

# Fishplot, 4 timepoints long sample
sample_list = ['PT_2_late', 'CTC_2_late', 'lung_2_late']
df_long = df_[sample_list].dropna()


df_long.sort_values('lung_2_late', ascending=False)


colors = { k:v for k,v in zip(df_long.index, sc.pl.palettes.godsnot_102[::-1])};
x = np.arange(df_long.shape[1])
y = [ np.array(x) for x in df_long.values.tolist() ]
grid = np.linspace(-1.3, 4, num=500)
y_smoothed = [ gaussian_smooth(x, y_, grid, .35) for y_ in y ]

# Viz
fig, ax = plt.subplots(figsize=(10,5))
ax.stackplot(grid, y_smoothed, baseline="sym", colors=colors.values())
format_ax(
    ax, xticks=sample_list, title='Longitudinal clones (n=75, late CTCs collection)',
    ylabel='Prevalence'
)
ax.spines[['right', 'bottom', 'top']].set_visible(False)
# ax.set_yticks([])
for l in x:    
    ax.axvline(x=l, color='k', linewidth=.5, linestyle='dashed')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'CTCs_long_dynamics.png'), dpi=400)


##
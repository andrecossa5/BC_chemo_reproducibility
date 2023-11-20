"""
Packed circles.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from plotting_utils._utils import *
from plotting_utils._colors import *
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
from BC_chemo_utils.plotting import *
matplotlib.use('macOSX')


#


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal')

# Read cells meta
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
meta = adata.obs[['GBC', 'sample', 'origin', 'condition', 'dataset']]

# Packed circles clones single-cell
df_freq = (
    meta.groupby(['condition', 'origin', 'sample'])
    ['GBC'].value_counts(normalize=True).loc[lambda x:x>0]
    .reset_index(name='freq')
    .rename(columns={'level_3':'GBC'})
    .merge(
        meta[['sample', 'dataset']]
        .reset_index(drop=True)
        .drop_duplicates(),
        
    )
    .drop_duplicates()
)

# Colors
clones_colors = (
    pd.read_csv(
        os.path.join(path_data, 'clones_colors_sc.csv'), 
        index_col=0
    )['color']
    .to_dict()
)


##


# Fig 
fig = plt.figure(figsize=(12.5, 5))

order = [
    'NT_NT_1', 'AC_NT_1', 'NT_AC_1', 'AC_AC_1', 'NT_NT_2', 'AC_NT_2',
    'NT_AC_2', 'AC_AC_2', 'NT_NT_3', 'AC_NT_3', 'NT_AC_3', 'AC_AC_3'
]

for i, dataset in enumerate(order):

    ax = fig.add_subplot(3, 4, i+1)
    axins1 = inset_axes(ax, width="40%", height="98%", loc='upper left')
    axins2 = inset_axes(ax, width="40%", height="98%", loc='upper right')
    df_ = df_freq.query('dataset==@dataset and freq>.01 and origin=="PT"').set_index('GBC')
    packed_circle_plot(
        df_, covariate='freq', ax=axins1, color=clones_colors, cmap='Reds', annotate=False
    )
    df_ = df_freq.query('dataset==@dataset and freq>.01 and origin=="lung"').set_index('GBC')
    packed_circle_plot(
        df_, covariate='freq', ax=axins2, color=clones_colors, cmap='Reds', annotate=False
    )
    ax.set(title=dataset)
    fig.tight_layout()
    ax.axis('off')

fig.subplots_adjust(bottom=.1, top=.9, left=.1, right=.9, wspace=.1, hspace=.5)
fig.savefig(os.path.join(path_results, 'circle_packed.png'), dpi=500)


##
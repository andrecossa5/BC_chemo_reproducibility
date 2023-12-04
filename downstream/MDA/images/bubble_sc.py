"""
Script to produce pro-metastatic clones summary and visualizations
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import random
from plotting_utils._utils import *
from plotting_utils._colors import *
from plotting_utils._plotting_base import *
from plotting_utils._plotting import *
matplotlib.use('macOSX')


#


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'clonal')

# Read cells meta
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
df = adata.obs[['GBC', 'sample', 'origin', 'condition']]

# Bubble plot single-cell
df_freq = (
    df.groupby(['condition', 'origin', 'sample'])
    ['GBC'].value_counts(normalize=True).loc[lambda x:x>0]
    .reset_index(name='freq')
    .rename(columns={'level_3':'GBC'})
)

# Random colors for clones
df_freq['sample'] = pd.Categorical(
    df_freq['sample'], 
    categories=df_freq['sample'].cat.categories[::-1]
)
clones = df_freq['GBC'].cat.categories

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

clones_colors = (
    pd.read_csv(
        os.path.join(path_data, 'clones_colors_sc.csv'), 
        index_col=0
    )['color']
    .to_dict()
)
df_freq['area_plot'] = df_freq['freq'] * (3000-5) + 5

# Fig 
fig, ax = plt.subplots(figsize=(6, 6))
scatter(df_freq, 'GBC', 'sample', by='GBC', c=clones_colors, s='area_plot', a=0.5, ax=ax)
format_ax(ax, title='Clones by sample', xlabel='Clones', xticks='')
ax.text(.4, .20, f'n cells total: {df.shape[0]}', transform=ax.transAxes)
ax.text(.4, .17, f'n clones total: {df_freq["GBC"].cat.categories.size}', transform=ax.transAxes)
s = df_freq.groupby('sample').size()
ax.text(.4, .14, f'n clones per sample: {round(s.median(),2)} (+-{round(s.std(),2)})', transform=ax.transAxes)
s = df_freq.query('origin == "PT"').groupby('sample')['freq'].median()
s = s[~s.isna()]
ax.text(.4, .11, f'Prevalence PT: {round(s.median(),2)} (+-{round(s.std(),2)})', transform=ax.transAxes)
s = df_freq.query('origin == "lung"').groupby('sample')['freq'].median()
s = s[~s.isna()]
ax.text(.4, .08, f'Prevalence mets: {round(s.median(),2)} (+-{round(s.std(),2)})', transform=ax.transAxes)
fig.tight_layout()


fig.savefig(os.path.join(path_results, 'bubbles_sc.png'), dpi=300)


##


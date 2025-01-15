"""
Script to produce pro-metastatic clones summary and visualizations
"""

import os
import pickle
import pandas as pd
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
df = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)

# Calculate clonal prevalences
df_freq = (
    df[['GBC', 'sample']]
    .groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .reset_index(name='freq')
    .rename(columns={'level_3':'GBC'})
    .merge(df[['sample', 'origin']].drop_duplicates(), on='sample')    
)

# Reorder samples
categories = [
    'NT_NT_PTs_1', 'NT_NT_mets_1', 'NT_NT_PTs_2', 'NT_NT_mets_2', 'NT_NT_PTs_3', 'NT_NT_mets_3',
    'NT_AC_PTs_1', 'NT_AC_mets_1', 'NT_AC_PTs_2', 'NT_AC_mets_2', 'NT_AC_PTs_3', 'NT_AC_mets_3',
    'AC_NT_PTs_1', 'AC_NT_PTs_2', 'AC_NT_PTs_3',
    'AC_AC_PTs_1', 'AC_AC_mets_1', 'AC_AC_PTs_2', 'AC_AC_mets_2', 'AC_AC_PTs_3', 'AC_AC_mets_3'
][::-1]
df_freq['sample'] = pd.Categorical(df_freq['sample'], categories=categories)
df_freq.sort_values(by=['sample'], inplace=True)

# Random colors for clones
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

# Read colors
with open(os.path.join(path_data, 'clones_colors_sc.pickle'), 'rb') as f:
    clones_colors = pickle.load(f)

# Set clones size
df_freq['area_plot'] = df_freq['freq'] * (3000-5) + 5

# Fig 
fig, ax = plt.subplots(figsize=(6, 6))
scatter(df_freq, 'GBC', 'sample', by='GBC', c=clones_colors, s='area_plot', a=0.5, ax=ax)
format_ax(ax, title='Clones by sample', xlabel='Clones', xticks='')
ax.text(.4, .20, f'n cells total: {df.shape[0]}', transform=ax.transAxes)
ax.text(.4, .17, f'n clones total: {df_freq["GBC"].unique().size}', transform=ax.transAxes)
s = df_freq.groupby('sample').size()
ax.text(.4, .14, f'n clones per sample: {round(s.median(),2)} (+-{round(s.std(),2)})', transform=ax.transAxes)
s = df_freq.query('origin == "PT"').groupby('sample')['freq'].median()
s = s[~s.isna()]
ax.text(.4, .11, f'Prevalence PT: {round(s.median(),2)} (+-{round(s.std(),2)})', transform=ax.transAxes)
s = df_freq.query('origin == "lung"').groupby('sample')['freq'].median()
s = s[~s.isna()]
ax.text(.4, .08, f'Prevalence mets: {round(s.median(),2)} (+-{round(s.std(),2)})', transform=ax.transAxes)
fig.tight_layout()

# Save
fig.savefig(os.path.join(path_results, 'bubbles_sc.png'), dpi=300)


##


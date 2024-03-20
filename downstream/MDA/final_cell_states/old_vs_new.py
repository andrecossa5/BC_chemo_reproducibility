"""
Relationship old cell_states vs new top clustering solution.
"""

import os
import numpy as np
import pandas as pd
from plotting_utils._plotting_base import *
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
matplotlib.use('macOSX')


##


# Set paths
path_data = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/MDA/cell_states/'

# Read info
old_anno = pd.read_csv(os.path.join(path_data, 'clusters_annotation_definitive.csv'), index_col=0)
df_old = pd.read_csv(os.path.join(path_data, 'old_clusters.csv'), index_col=0)
df_new = pd.read_csv(os.path.join(path_data, 'new_clusters.csv'), index_col=0)

# Relationship
common = set(df_old.index) & set(df_new.index)

# Go
df = (
    df_old.loc[list(common)].reset_index()
    .merge(old_anno.reset_index(), on='50_NN_0.94')
    .set_index('index')
    .join(df_new.loc[list(common)][['leiden']])
)

# Concordance
normalized_mutual_info_score(df['Annotation_definitive'], df['leiden'])
adjusted_rand_score(df['Annotation_definitive'], df['leiden'])

# Correspondence in a cont table
cont = pd.crosstab(df['Annotation_definitive'], df['leiden'], normalize=1)
cont.to_csv(os.path.join(path_data, 'contingency_table.csv'))


##

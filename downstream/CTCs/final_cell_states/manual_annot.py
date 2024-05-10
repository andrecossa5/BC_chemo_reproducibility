"""
Manual annotation of cell states.
"""

import os
import pickle
from Cellula.plotting._colors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from BC_chemo_utils.tests import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'CTCs')


##


# Data
# adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
embs = pd.read_csv(os.path.join(path_data, 'full_embs.csv'), index_col=0)
cell_state_map = (
        pd.read_csv(
        os.path.join(path_data, 'CTCs_cluster_annotation.csv'), 
        index_col=0
    )['Cell state'].to_dict()
)


##


# New, final annot
# adata.obs['leiden'] = adata.obs['leiden'].astype('int')
# adata.obs['final_cell_state'] = adata.obs['leiden'].map(cell_state_map)

# Save clustered and full_embs
# adata.write(os.path.join(path_data, 'clustered.h5ad'))
embs['final_cell_state'] = embs['leiden'].map(cell_state_map)
embs.to_csv(os.path.join(path_data, 'full_embs.csv'))


##
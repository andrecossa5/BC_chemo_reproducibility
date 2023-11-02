"""
CosPar.
"""

import os
import re
import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from anndata import AnnData
import cospar as cs
from Cellula.preprocessing._pp_recipes import *
from Cellula.preprocessing._pp import *
from Cellula.preprocessing._embeddings import *
from Cellula.preprocessing._neighbors import *
from Cellula.preprocessing._integration import *
from Cellula.clustering._clustering import *
from Cellula.plotting._colors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results', 'embeddings')

##

# Read 
adata = sc.read(os.path.join(path_data, 'subsampled_with_embs.h5ad'))
# Format
a = cs.pp.initialize_adata_object(
    X_state=adata.layers['raw'],
    gene_names=adata.var_names,
    cell_names=adata.obs_names, 
    time_info=adata.obs['origin'],
    state_info=adata.obs['leiden'],
    X_pca=adata.obsm['X_harmony'],
    X_emb=adata.obsm['X_umap'],
    data_des='MDA_chemo'
)
df_clone = adata.obs.loc[:, ['GBC']].reset_index().rename(columns={'index':'Cell_ID', 'GBC':'Clone_ID'})
cs.pp.get_X_clone(a, df_clone["Cell_ID"], df_clone["Clone_ID"].astype('str'))


##


################################## Clonal analysis
selected_times = None
cs.pl.barcode_heatmap(
    a,
    color_bar=True,
    log_transform=False,
    binarize=True,
    plot=True
)

# Fate coupling
cs.tl.fate_coupling(a, source='X_clone')  # compute the fate coupling
cs.pl.fate_coupling(a, source='X_clone') 

# Fate hierarchy
cs.tl.fate_hierarchy(a, source='X_clone')  # compute the fate hierarchy
cs.pl.fate_hierarchy(a, source='X_clone')

# Clonal fate bias
cs.tl.clonal_fate_bias(
    a, selected_fate=[2], alternative="two-sided"
)  # compute the fate hierarchy
cs.pl.clonal_fate_bias(a) 

# Inspect
ids = a.uns['clonal_fate_bias']["Clone_ID"][:2]
cs.pl.clones_on_manifold(
    a,
    selected_clone_list=ids,
    color_list=["black", "red", "blue"],
    clone_markersize=2,
)

##################################


##


################################## Transition map inference
a = cs.tmap.infer_Tmap_from_multitime_clones(
    a,
    clonal_time_points=["PT"],
    later_time_point="lung",
    smooth_array=[20, 15, 10, 5],
    sparsity_threshold=0.1,
    intraclone_threshold=0.2,
    max_iter_N=10,
    epsilon_converge=0.01,
)
# a.write(os.path.join(path_main, 'results', 'coSpar', 'adata_with_tmap.h5ad'))

##

# Visualization of single-cell phate maps
cs.pl.single_cell_transition(
    a,
    selected_state_id_list=[2],
    source="intraclone_transition_map",
    map_backward=True,
    savefig=True
)

# Visualization of states
cs.tl.fate_map(
    a,
    selected_fates=[0,1],
    source="intraclone_transition_map",
    map_backward=True,
)
cs.pl.fate_map(
    a,
    selected_fates=[1],
    source="intraclone_transition_map",
    plot_target_state=True
)
plt.show()

# Fate potency
cs.tl.fate_potency(
    a,
    source="intraclone_transition_map",
    map_backward=True,
    method="norm-sum",
    fate_count=True,
)
cs.pl.fate_potency(a, source="intraclone_transition_map")

##

"""
Perturbation analysis.
"""

import os
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import scgen
import pertpy as pt
import matplotlib
from plotting_utils import *

matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'perturbation')

##

# Read clustered adata
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))


##


# Augur
ag_rfc = pt.tl.Augur("random_forest_classifier")
loaded_data = ag_rfc.load(adata, label_col="condition", cell_type_col="leiden")
v_adata, v_results = ag_rfc.predict(
    loaded_data, subsample_size=20, n_threads=8, select_variance_features=True, span=1
)
v_results["summary_metrics"]


# scgen
red = adata[:, adata.var['highly_variable_features']].copy()
scgen.SCGEN.setup_anndata(red, batch_key="condition", labels_key="leiden")
model = scgen.SCGEN(red, n_hidden=800, n_latent=100, n_layers=2)
model.train(
    max_epochs=100, batch_size=32, early_stopping=True, early_stopping_patience=25
)


# Visualize trends



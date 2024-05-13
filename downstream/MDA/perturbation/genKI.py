"""
Perturbation analysis.
"""

import os
import scanpy as sc
import pandas as pd
import GenKI as gk
from scipy.sparse import csr_matrix
from GenKI.dataLoader import DataLoader
from GenKI.train import VGAE_trainer
from GenKI import utils
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'perturbation')

##

# Read clustered adata
adata = sc.read(os.path.join(path_data, 'subsample.h5ad'))

# Select HVGs, and scale 
genes = adata.var['highly_variable_features']
adata = adata[:,genes].copy()
adata.layers['norm'] = csr_matrix(adata.X)
sc.pp.scale(adata)
adata.X = csr_matrix(adata.X)
# adata = build_adata('/Users/IEO5505/Downloads/microglial_seurat_WT.h5ad')

# GenKI

# Load
gene = 'PAEP'
data_wrapper =  DataLoader(adata, target_gene=[gene], rebuild_GRN=False, pcNet_name="pcNet", n_cpus=8)
data_wt = data_wrapper.load_data()
data_ko = data_wrapper.load_kodata()

# Train
hyperparams = {"epochs": 100, "lr": 7e-4,  "beta": 1e-4, "seed": 1234}
sensei = VGAE_trainer(
    data_wt, epochs=hyperparams["epochs"], lr=hyperparams["lr"], 
    log_dir=None, beta=hyperparams["beta"], seed=hyperparams["seed"]
)
sensei.train()

# Rank KO-affected genes
z_mu_wt, z_std_wt = sensei.get_latent_vars(data_wt)
z_mu_ko, z_std_ko = sensei.get_latent_vars(data_ko)
dis = gk.utils.get_distance(z_mu_ko, z_std_ko, z_mu_wt, z_std_wt, by="KL")
res_raw = utils.get_generank(data_wt, dis, rank=True)
res_raw

# Get significantly affected by permutation test
null = sensei.pmt(data_ko, n=100, by="KL")
res = utils.get_generank(data_wt, dis, null)
res.head(100).to_csv('prova.csv')




##
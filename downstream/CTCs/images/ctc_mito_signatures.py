"""
Script to analyze mitochondrial dynamics and mechanotransduction signatures
"""

import os
import base64
import zlib
import json
import umap
import matplotlib.axes
import textalloc as ta
from typing import Dict, Iterable, Any, Tuple
from matplotlib.lines import Line2D 
import anndata as ad
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns
import pickle 
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.stats import spearmanr
from scipy.sparse import csr_matrix
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *
from Cellula.dist_features._Gene_set import *
from Cellula.dist_features._signatures import scanpy_score, wot_zscore
from Cellula.preprocessing._neighbors import *
from Cellula.clustering._clustering import *
from Cellula.plotting._plotting import *
from matplotlib.gridspec import GridSpec
from plotting_utils import plotting_base
import plotting_utils as plu
from BC_chemo_utils.utils import prep_df_for_dotplot
from BC_chemo_utils.tests import *
matplotlib.use('macOSX')



#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility"
path_data= os.path.join(path_main,"data", "CTCs")
path_results= os.path.join(path_main, "results", "CTCs","mito_signatures")

#Data
adata= sc.read_h5ad(os.path.join(path_data, "clustered.h5ad"))
sign_fission = pd.read_csv(os.path.join(path_data, "GOBP_MITOCHONDRIAL_FISSION.v2025.1.Hs.tsv"), sep='\t')
sign_fusion = pd.read_csv(os.path.join(path_data, "GOBP_MITOCHONDRIAL_FUSION.v2025.1.Hs.tsv"), sep='\t')
sign_mecha = ['YAP1','WWTR1','TAZ','MFN1','MFN2','OPA1','FIS1','PRKCE','YME1L1','LETM1','RHOT1','RHOT2','PGC1A',
              'PGC1B','TFAM','SRC','TOMM20','RHOA','ITGB1','TOMM40','HSPD1', 'ROCK1','ROCK2']
sign_oxphos = pd.read_csv(os.path.join(path_data, "GOBP_OXIDATIVE_PHOSPHORYLATION.v2025.1.Hs.tsv"), sep='\t')

fission = sign_fission['GOBP_MITOCHONDRIAL_FISSION'][16].split(',')
fusion = sign_fusion['GOBP_MITOCHONDRIAL_FUSION'][16].split(',')
oxphos = sign_oxphos['GOBP_OXIDATIVE_PHOSPHORYLATION'][16].split(',')
mito_dynamics = list(set(fission + fusion))

dynamics_present = []
for gene in mito_dynamics:
    if gene in adata.var_names:
        dynamics_present.append(gene)

mecha_present = []
for gene in sign_mecha:
    if gene in adata.var_names:
        mecha_present.append(gene)

oxphos_present = []
for gene in oxphos:
    if gene in adata.var_names:
        oxphos_present.append(gene)

#For Cellula
with open(os.path.join(path_data, 'dynamics_present.txt'), 'w') as f:
    for gene in dynamics_present:
        f.write(f"dynamics_present\t{gene.strip()}\n")

with open(os.path.join(path_data, 'mecha_present.txt'), 'w') as f:
    for gene in mecha_present:
        f.write(f"mecha_present\t{gene.strip()}\n")

with open(os.path.join(path_data, 'oxphos_present.txt'), 'w') as f:
    for gene in oxphos_present:
        f.write(f"oxphos_present\t{gene.strip()}\n")


with open(os.path.join(path_data,'signatures.pickle'), 'rb') as f:
    data = pickle.load(f)

#boxplot lung CTCs signatures
signatures = data['scores']
adata.obs['dynamics_present'] = signatures['dynamics_present']
adata.obs['mecha_present'] = signatures['mecha_present']
df = pd.DataFrame({
    'dynamics_present' : adata.obs['dynamics_present'],
    'mecha_present' : adata.obs['mecha_present'],
    'origin': adata.obs['origin']
})
df_subset = df[df['origin'].isin(['lung', 'ctc'])].copy()

fig, ax = plt.subplots(figsize=(6.9,5))
order = ['lung','CTC']

plu.box(df_subset,
    x='origin', y='dynamics_present', ax=ax, color='grey', add_stats=True, 
    pairs=[['lung', 'CTC']], 
    x_order=order
)
plu.strip(df_subset, x='origin', y='dynamics_present', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mitochondrial dynamic signature in lung and ctc', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_dynamics_signature.png'), dpi=300)

#scanpy_score
gene_signatures = {
    'mito_dynamics': dynamics_present,
    'mechano_mito': mecha_present,
    'oxphos': oxphos_present
}

import numpy as np
signature_names = list(gene_signatures.keys())

for i, name in enumerate(signature_names):
    genes = gene_signatures[name]
    scores = scanpy_score(adata, genes)
    adata.obs[f'{name}_score'] = scores.values.flatten() 

df = pd.DataFrame({
    'sample': adata.obs['sample'],
    'dynamics_present' : adata.obs['mito_dynamics_score'],
    'mecha_present' : adata.obs['mechano_mito_score'],
    'oxphos': adata.obs['oxphos_score'],
    'origin': adata.obs['origin']
})

df_subset = df[df['origin'].isin(['lung', 'CTC'])].copy()
df_pseudo = df_subset.groupby(['sample', 'origin'])[['mecha_present']].mean().reset_index()


fig, ax = plt.subplots(figsize=(6.9,5))
order = ['lung','CTC']

plu.box(df_valid,
    x='origin', y='mecha_present', ax=ax, color='grey', add_stats=True, 
    pairs=[['lung', 'CTC']], 
    x_order=order
)
plu.strip(df_pseudo, x='origin', y='mecha_present', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mitochondrial mechanotransduction signature in lung and ctc', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_mecha_signature_pseudobulk.png'), dpi=300)


#UMAP of oxphos coupled with origin
fig, axs = plt.subplots(1, 2, figsize=(10, 5))

# oxphos score with vmax
sc.pl.umap(
    adata,
    color='oxphos_score',
    cmap='viridis',
    vmax=1,
    size=5,
    ax=axs[0],
    title='Oxidative phosphorylation signature',
    show=False,
)

sc.pl.umap(
    adata,
    color='origin',
    cmap='viridis',
    vmax=2,
    size=5,
    ax=axs[1],
    show=False,
    title='Origin'
)

plt.tight_layout()
plt.savefig(os.path.join(path_results, 'umap_oxphos_origin'), dpi=300)

#UMAP for pro-CTC 
pro_ctc_barcodes = [
    "TATGGCGGTTGTGCTGTC", "CCATAAGGGGAGGATGTA", "TGACTGTGGAGTCCTGAA",
    "TTTGTACGTGGGACCATA", "CGCGACTAGAGACGGCTC", "CTGCGGTTTCGTTAACGC",
    "ATCTGGGACCTCAACAAG", "CGCGTCACACTGTCGGGC", "CTCGAGTCTGGATGCGTA",
    "GGTGCAGGCGGTGTGGGG", "TTCCGAGCAGCTGCTTGG"
]
adata.obs["pro"] = "Other"
is_pro_ctc = adata.obs["GBC"].isin(pro_ctc_barcodes)
adata.obs.loc[is_pro_ctc, "pro"] = "pro_ctc_" + adata.obs.loc[is_pro_ctc, "origin"].astype(str)
adata.obs["oxphos_pro_only"] = np.nan
adata.obs.loc[adata.obs["pro"] != "Other", "oxphos_pro_only"] = adata.obs.loc[adata.obs["pro"] != "Other", "oxphos_score"]

adata.obs["origin_pro_only"] = np.nan
adata.obs.loc[adata.obs["pro"] != "Other", "origin_pro_only"] = adata.obs.loc[adata.obs["pro"] != "Other", "origin"].astype(str)

# Set proper order
adata.obs["origin_pro_only"] = pd.Categorical(
    adata.obs["origin_pro_only"],
    categories=["PT", "lung", "CTC"]
)

fig,ax = plt.subplots(figsize=(5, 5))

sc.pl.umap(
    adata,
    color='oxphos_pro_only',
    cmap='viridis',
    vmax=1,
    size=5,
    ax=ax,
    title='OXPHOS (pro-CTC clones only)',
    show=False
)

plt.tight_layout()
plt.savefig(os.path.join(path_results, 'umap_oxphos_pro_ctc'), dpi=300)

fig,ax = plt.subplots(figsize=(5.5, 5))

sc.pl.umap(
    adata,
    color='origin_pro_only',
    cmap='viridis',
    vmax=1,
    size=5,
    ax=ax,
    title='Origin (pro-CTC clones only)',
    show=False
)

plt.tight_layout()
plt.savefig(os.path.join(path_results, 'umap_origin_pro_ctc'), dpi=300)



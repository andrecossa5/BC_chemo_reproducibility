"""
Script to analyze pyscenic results CTC
"""

import os
import base64
import zlib
import json
import umap
import matplotlib.axes
from gseapy import ssgsea
import textalloc as ta
from typing import Dict, Iterable, Any, Tuple
from matplotlib.lines import Line2D 
from scipy.cluster.hierarchy import leaves_list, linkage
import anndata as ad
import pandas as pd
import scanpy as sc
import loompy as lp
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse import csr_matrix
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *
from Cellula.dist_features._Gene_set import *
from Cellula.dist_features._signatures import scanpy_score, wot_zscore
from Cellula.preprocessing._neighbors import *
from Cellula.clustering._clustering import *
from matplotlib.gridspec import GridSpec
from plotting_utils import plotting_base
import plotting_utils as plu
from BC_chemo_utils.utils import prep_df_for_dotplot
from BC_chemo_utils.tests import *
matplotlib.use('macOSX')




#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility"
path_data= os.path.join(path_main,"data", "CTCs","grn","resources")
path_results= os.path.join(path_main, "results", "CTCs", "gene_reg_net")

#Data
df_motif= pd.read_csv(os.path.join(path_main, "data","CTCs", "grn","resources", "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),sep="\t")
regulons= pd.read_csv(os.path.join(path_data, "ctc_regulons.csv"))
adata= sc.read_h5ad(os.path.join(path_data, "clustered.h5ad"))
adata_sc = sc.read_h5ad(os.path.join(path_data, "clustered_norm.h5ad"))
lf= lp.connect(os.path.join(path_data, "ctc_aucell.loom"), mode="r+", validate=False)
auc_mtx= pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
rank_df_10 = pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
rank_df_500 = pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
score_df_10=pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"))
score_df_500= pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"))
#reformat motif 
df_motif.columns = df_motif.columns.str.lstrip('#*')

#extract metadata 
encoded= lf.attrs['MetaData']
compressed_data= base64.b64decode(encoded)
decompressed_data = zlib.decompress(compressed_data)
metadata = json.loads(decompressed_data)

#Reformat regulons 
regulons.columns = regulons.iloc[0]         
regulons = regulons.drop(index=[0])     
regulons.columns.values[0] = regulons.iloc[0, 0]  
regulons.columns.values[1] = regulons.iloc[0, 1]
regulons= regulons.drop(index=[1]).reset_index(drop=True)


#New Anndata
common_cells = auc_mtx.index.intersection(adata.obs_names)
obs_subset = adata.obs.loc[common_cells].copy()
auc_mtx = auc_mtx.loc[common_cells]

adata_auc = ad.AnnData(X=auc_mtx.values, obs=obs_subset, var=pd.DataFrame(index=auc_mtx.columns))
auc_mtx.columns = auc_mtx.columns.str.replace('(+)', '', regex=False)
adata_auc.var['regulon'] = auc_mtx.columns
adata_auc.var_names = adata_auc.var['regulon']

def calculate_percent_cells_for_regulons(adata_auc):
    """
    Calculate the percentage of cells expressing each regulon.
    Assumes that regulons are in `adata.var` and the expression data is in `adata.X`.
    """
    regulon_expr_matrix = adata_auc.X > 0
    percent_cells = np.sum(regulon_expr_matrix, axis=0) / regulon_expr_matrix.shape[0] * 100

    adata_auc.var['percent_cells'] = percent_cells

calculate_percent_cells_for_regulons(adata_auc)

#add highly_variable_features
highly_variable_features = adata.var['highly_variable_features']
regulon_to_gene_mapping = {}  
for regulon_name in adata_auc.var.index:
    gene_name = regulon_name  
    regulon_to_gene_mapping[regulon_name] = gene_name

adata_auc.var['highly_variable_features'] = [
    highly_variable_features.get(regulon_to_gene_mapping[regulon], False)  
    for regulon in adata_auc.var.index
]

# add layer RAW 
adata_auc.layers['raw'] = adata_auc.X

#add 'mean' and 'var' columns
mean_values = np.array(adata_auc.X.mean(axis=0)).flatten()
variance_values = np.array(adata_auc.X.var(axis=0))

adata_auc.var['mean'] = mean_values
adata_auc.var['var'] = variance_values

X = adata_auc.X

# Check if the matrix is dense (numpy ndarray), and if so, convert it to sparse CSR matrix
if isinstance(X, np.ndarray):
    print("Converting dense matrix to sparse CSR format...")
    adata_auc.X = csr_matrix(X)

if adata_auc.X.dtype != 'float32':
    adata_auc.X = adata_auc.X.astype('float32')

#distribution
if sparse.issparse(adata_auc.X):
    data = adata_auc.X.toarray().flatten()
else:
    data = adata_auc.X.flatten()
plt.figure(figsize=(8,5))
plt.hist(data, bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("AUCell scores distribution")
plt.tight_layout()
plt.savefig(os.path.join(path_results, "distribution_aucell.png"),dpi=300)


#Regulon : tf + gene_set
alt_mtx = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)


#Create dictionary 
d_reg = {}

if isinstance(metadata.get('regulonThresholds', None), list):
    for item in metadata['regulonThresholds']:
        tf = item.get('regulon')
        motif_data = item.get('motifData', {})
        genes_for_set = alt_mtx.index[alt_mtx[tf] == 1].tolist()
        d_reg[tf] = {
            'gene_set': genes_for_set,
            'motif_data': motif_data
        }

#save d_reg
rows = []
for tf, data in d_reg.items():
    gene_set = data['gene_set']
    motif_data = data['motif_data']
    
    rows.append({
        'regulon': tf,
        'gene_set': ', '.join(gene_set),  
        'motif_data': str(motif_data)    
    })

df = pd.DataFrame(rows)
df.to_csv(os.path.join(path_results,'ctc_regulon_data.csv'), index=False)

#re-score values 
regulon_name = list(d_reg.keys())
n_cells= adata.n_obs
X = np.zeros((n_cells,len(d_reg)))
for i , key in enumerate(regulon_name):
    regulon=d_reg[key]['gene_set']
    scores = scanpy_score(adata,regulon)
    X[:,i] = scores.values.flatten()

#create matrix 
X_sparse= csr_matrix(X)

#create new Anndata
common_cells = auc_mtx.index.intersection(adata.obs_names)
obs_subset = adata.obs.loc[common_cells].copy()
auc_mtx = auc_mtx.loc[common_cells]

adata_sc = ad.AnnData(X=X_sparse, obs=obs_subset, var=pd.DataFrame(index=regulon_name))

X_dense= adata_sc.X.toarray()

plt.figure(figsize=(8,5))
plt.hist(X_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("Distribution of All Scores")
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(path_results, "distribution_scanpyscores.png"), dpi=300)

#normalized scanpy score 
Z = adata_sc.X

if sparse.issparse(Z):
    Z = Z.toarray()

X_zscored= (Z - Z.mean(axis=0))/ Z.std(axis=0)

adata_sc.X = csr_matrix(X_zscored)
adata_sc.write_h5ad(os.path.join(path_data,'clustered_norm.h5ad'))
X_dense= adata_sc.X.toarray()

plt.figure(figsize=(8,5))
plt.hist(X_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("scanpy distribution (Z-score normalized)")
plt.tight_layout()
plt.savefig(os.path.join(path_results, "distribution_scanpy_zscored.png"), dpi=300)

#add necessary columns for DE 
#create column percent_cells
regulon_expr_matrix = adata_sc.X > 0
percent_cells = np.sum(regulon_expr_matrix, axis=0) / regulon_expr_matrix.shape[0] * 100
adata_sc.var['percent_cells'] = np.ravel(percent_cells)

# add layer RAW 
adata_sc.layers['raw'] = adata_sc.X

#add 'mean' and 'var' columns
T = adata_sc.X

if sparse.issparse(T):
    T = T.toarray()

adata_sc.X = T
mean_values = np.array(adata_sc.X.mean(axis=0))
variance_values = np.array(adata_sc.X.var(axis=0))

adata_sc.var['mean'] = mean_values
adata_sc.var['var'] = variance_values

#regulon in var
auc_mtx.columns = auc_mtx.columns.str.replace('(+)', '', regex=False)
adata_sc.var['regulon'] = auc_mtx.columns
adata_sc.var_names = adata_sc.var['regulon']
#put all regulons as highly_variable_features
adata_sc.var['highly_variable_features'] = adata_sc.var['regulon']

# Check if the matrix is dense, if so, convert it to sparse CSR matrix
if isinstance(T, np.ndarray):
    print("Converting dense matrix to sparse CSR format...")
    adata_sc.X = csr_matrix(T)

if adata_sc.X.dtype != 'float32':
    adata_sc.X = adata_sc.X.astype('float32')

## DE ##

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata_sc, path_data, contrasts_name='pt_ctc')

# Here we go
D = Dist_features(adata_sc, contrasts, jobs=jobs)
D.run_all_jobs()

dfs = []
categories= ['pt_ctc','lung_ctc','pt_lung']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, 'PT_CTC_DEGs.csv'))

promet=all_degs[all_degs['comparison']== 'promet_AC_vs_NT0_vs_promet_AC_vs_NT1']
promet['effect_size']
promet.to_csv(os.path.join(path_results,'DEgs_promet_AC_vs_promet_NT.csv'))
all_degs['comparison'].unique()


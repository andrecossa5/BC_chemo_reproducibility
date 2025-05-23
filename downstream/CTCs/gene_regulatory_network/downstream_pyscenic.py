"""
Script to analyze pyscenic results
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
path_data= os.path.join(path_main,"data", "MDA")
path_results= os.path.join(path_main, "results", "MDA")

#Data
df_motif= pd.read_csv(os.path.join(path_main, "data","CTCs", "grn","resources", "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),sep="\t")
regulons= pd.read_csv(os.path.join(path_data, "PAEP_regulons_dropout.csv"))
adata= sc.read_h5ad(os.path.join(path_data, "clustered.h5ad"))
lf= lp.connect(os.path.join(path_data, "PAEP_aucell_dropout.loom"), mode="r+", validate=False)
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

#new files (d_reg, adata_auc_thrb, adata_sc)
d_reg= pd.read_csv(os.path.join(path_results,'regulon_data.csv'))
adata_auc_thrb= sc.read(os.path.join(path_data,"clustered_thrb.h5ad"))
adata_sc = sc.read(os.path.join(path_data,"clustered_scanpy_no_norm.h5ad"))
# adata_sc.layers['raw'] = adata_sc.X
# adata_sc.var['regulon'] = adata_sc.var_names
# adata_sc.var['regulon'] = adata_sc.var_names.str.replace('(+)', '', regex=False)
# adata_sc.var_names = adata_sc.var['regulon']
# adata_sc.obs['THRB_score'] = adata_sc[:,"THRB"].X.toarray().flatten()
# adata_sc.obs['comparison'] = adata_auc_thrb.obs['comparison']
# adata_sc.write(os.path.join(path_data, "clustered_norm.h5ad"))
adata_sc=sc.read(os.path.join(path_data, "clustered_norm.h5ad"))
adata_sc.obs.drop('THRB_score', axis=1)
all_degs=pd.read_csv(os.path.join(path_data, "Degs_regulon_scanpy_score.csv"))


#create new Anndata
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
plt.show()


#stats ptx
auc_mtx
gene_counts = [len(d_reg[key]['gene_set']) for key in d_reg]
mean_gene_count = sum(gene_counts) / len(gene_counts)
min_gene = min(gene_counts)
max_gene = max(gene_counts)




#normalized aucell score 
Z = adata_auc.X

if sparse.issparse(Z):
    Z = Z.toarray()

X_zscored= (Z - Z.mean(axis=0))/ Z.std(axis=0)

adata_auc.Z = csr_matrix(X_zscored)

X_dense= adata_auc.Z.toarray()

plt.figure(figsize=(8,5))
plt.hist(X_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("AUCell distribution (Z-score normalized)")
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(path_results, "distribution_aucell_zscored.png"), dpi=3)







## DE ## AUCELL SCORES

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata_auc, path_data, contrasts_name='paep_contrasts')

# Here we go
D = Dist_features(adata_auc, contrasts, jobs=jobs)
D.select_genes()
for k in D.jobs:
    for x in D.jobs[k]:
        if x['model'] == 'wilcoxon':
            job_key = '|'.join([k, x['features'], x['model']])
            de_results, gene_set_dict = D.compute_DE(contrast_key=k, which='perc_1_no_miribo')
            D.Results.add_job_results(de_results, gene_set_dict, job_key=job_key)

dfs = []
categories= ['nonpro_AC_vs_NT', 'promet_AC_vs_NT', 'promet_AC', 'pro_nonpro_AC', 'promet_NT']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, "Degs_regulon.csv"))
#where THRB is DE 
degs_promet_AC=D.Results.results['promet_AC|genes|wilcoxon']['df']

  
df_markers = pd.read_csv(os.path.join(path_data, 'Degs_regulon.csv'), index_col=0)
df_markers['comparison'].unique()




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
df.to_csv(os.path.join(path_results,'regulon_data.csv'), index=False)


#Contrasts tested regulons
d_contrast={}
reg_to_extract= ['THRB(+)','SMAD5(+)','NR2F1(+)']
for key,value in d_reg.items():
    if key in  reg_to_extract:
        d_contrast[key] = value

#N regulons with PAEP
d_reg['gene_set'] = d_reg['gene_set'].apply(lambda x: [gene.strip() for gene in x.split(',')])
reg_paep_df = d_reg[d_reg['gene_set'].apply(lambda genes: 'PAEP' in genes)]
reg_paep = {
    row['regulon']: {
        'gene_set': row['gene_set'],
        'motif_data': row.get('motif_data', {})
    }
    for _, row in reg_paep_df.iterrows()
}




#Viz 
#create column comparison
contrasts
contrast= D.contrasts['pro_nonpro_AC']
adata_tmp.obs['pro_nonpro_AC'] = contrast.category
cols_to_merge = ['nonpro_AC_vs_NT', 'promet_AC_vs_NT', 'promet_AC', 'promet_NT']
adata_auc.obs['comparison'] = adata_auc.obs[cols_to_merge].apply(
    lambda row: next((val for val in row if val != 'to_exclude'), pd.NA),
    axis=1
)

print(adata_auc.obs[adata_auc.obs['comparison'].isna()])
adata_auc.obs['comparison'].unique()

adata_auc_thrb = adata_auc.copy()
thrb_val= auc_mtx['THRB']
adata_auc_thrb.obs['THRB'] = thrb_val

adata_auc_thrb.write(os.path.join(path_data, "clustered_thrb.h5ad"))
adata_auc_thrb= sc.read(os.path.join(path_data,"clustered_thrb.h5ad"))

#functions
def violin(
    df: pd.DataFrame, 
    x: str, 
    y: str, 
    by: str = None, 
    color: str = None,
    categorical_cmap: str|Dict[str,Any] = 'tab10', 
    x_order: Iterable[str] = None,
    add_stats: bool = False,
    pairs: Iterable[Iterable[str]] = None, 
    by_order: Iterable[str] = None,
    linewidth: float|str = .5, 
    ax: matplotlib.axes.Axes = None, 
    kwargs: Dict[str,Any] = {}
    ) -> matplotlib.axes.Axes:

    params = {   
        'inner' : 'quart'
    }    
    params = plu.update_params(params, kwargs)

    # Handle colors and by
    if by is None:
        sns.violinplot(
            data=df, x=x, y=y, ax=ax, 
            order=x_order, 
            color=color,
            linewidth=linewidth, 
            **params
        )
        
    elif by is not None and by in df.columns:
        # Categorical
        if pd.api.types.is_string_dtype(df[by]) or df[by].dtype == 'category':
            
            if isinstance(categorical_cmap, str):
                _cmap = plu.create_palette(df, by, palette=categorical_cmap)
            else:
                _cmap = categorical_cmap

            assert all([ x in _cmap for x in df[by].unique() ])
            sns.violinplot(
                data=df, x=x, y=y, ax=ax, 
                order=x_order, 
                dodge=True,
                hue=by, hue_order=by_order, palette=_cmap,
                linewidth=linewidth, 
                **params
            )
            ax.get_legend().remove()

        else:
            raise ValueError(f'{by} must be categorical or string!')
    
    else:
        raise KeyError(f'{by} not in df.columns!')

    if add_stats:
        plu.add_wilcox(df, x, y, pairs, ax, order=x_order)

    return ax

#to delete
adata_auc_thrb_clean = adata_auc_thrb[adata_auc_thrb.obs['comparison'].notna()].copy()

adata_auc_thrb_clean.obs['comparison'].unique()
#Here we go
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

pairs = [
    ['promet_NT', 'nonpro_NT'],
    ['promet_AC', 'nonpro_AC'],
    ['promet_AC', 'promet_NT'],
    ['nonpro_AC', 'nonpro_NT'],
    ['nonpro_AC','promet_NT']
]

violin(
    df=adata_auc_thrb_clean.obs,
    x='comparison',
    y='THRB',
    ax=ax,
    add_stats=True,
    pairs=pairs,
    linewidth=0.5
)

# Format the axis
plu.format_ax(ax, title='THRB scores', 
          xticks=adata_auc_thrb_clean.obs['comparison'].cat.categories, ylabel='Score', reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'THRB_violin.png'), dpi=400)


adata_auc_thrb.var.drop('THRB', inplace=True)

#umap
adata_mp=adata.copy()
adata_mp.obs['THRB_score'] = adata_sc[:,"THRB"].X.toarray().flatten()
adata_mp.obs['comparison']= adata_auc_thrb.obs['comparison']
regulon_to_plot=['THRB_score']
adata_mp.obs['promet_NT'] = adata_mp.obs['comparison'].apply(lambda x: x if x == 'promet_NT' else np.nan)
adata_mp.obs['THRB_score_promet_NT'] = adata_mp.obs.apply(
    lambda row: row['THRB_score'] if pd.notnull(row['promet_NT']) else np.nan,
    axis=1
)

#save umap
plt.figure(figsize=(6, 6))  
sc.pl.umap(adata_mp, color='comparison', cmap='viridis', vmin=0, show=False)
plt.savefig(os.path.join(path_results,"umap_comparison.png"), dpi=300, bbox_inches='tight')
plt.close()

plt.figure(figsize=(6, 6))  
sc.pl.umap(adata_mp, color=regulon_to_plot, cmap='viridis', vmin=0, show=False)
plt.show()
plt.savefig(os.path.join(path_results,"umap_THRB.png"), dpi=300, bbox_inches='tight')
plt.close()



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

adata_auc = ad.AnnData(X=X_sparse, obs=obs_subset, var=pd.DataFrame(index=regulon_name))

X_dense= adata_auc.X.toarray()

plt.figure(figsize=(8,5))
plt.hist(X_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("Distribution of All Scores")
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(path_results, "distribution_scanpyscores.png"), dpi=300)


# regulon_name = list(d_reg.keys())
# n_cells= adata.n_obs
# N = np.zeros((n_cells,len(d_reg)))
# for i , key in enumerate(regulon_name):
#     print(i)
#     regulon=d_reg[key]['gene_set']
#     scores = wot_zscore(adata,regulon)
#     N[:,i] = scores.values.flatten()

# #create matrix 
# N_sparse= csr_matrix(N)

# #create new Anndata
# common_cells = auc_mtx.index.intersection(adata.obs_names)
# obs_subset = adata.obs.loc[common_cells].copy()
# auc_mtx = auc_mtx.loc[common_cells]

# adata_reg = ad.AnnData(X=N_sparse, obs=obs_subset, var=pd.DataFrame(index=regulon_name))

# N_dense= adata_reg.X.toarray()

# plt.figure(figsize=(8,5))
# plt.hist(N_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
# plt.xlabel("Score")
# plt.ylabel("Frequency")
# plt.title("Distribution of All Scores")
# plt.tight_layout()
# plt.show()
# plt.savefig(os.path.join(path_results, "distribution_wot_zscore.png"), dpi=300)



#z-score normalization
W = adata_sc.X

if sparse.issparse(W):
    W = W.toarray()

X_zscored= (W - W.mean(axis=0))/ W.std(axis=0)

adata_sc.X = csr_matrix(X_zscored)

X_dense= adata_sc.X.toarray()

plt.figure(figsize=(8,5))
plt.hist(X_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("Scanpy scores distribution (Z-score normalized)")
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(path_results, "distribution_scanpy_scorezscored.png"), dpi=300)


# T = adata_reg.X

# if sparse.issparse(T):
#     T = T.toarray()

# T_zscored = (T - T.mean(axis=0))/ T.std(axis=0)
# adata_reg.X = csr_matrix(T_zscored)

# N_dense= adata_reg.X.toarray()

# plt.figure(figsize=(8,5))
# plt.hist(N_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
# plt.xlabel("Score")
# plt.ylabel("Frequency")
# plt.title("Distribution of All Scores")
# plt.tight_layout()

# plt.savefig(os.path.join(path_results, "distribution_wotzscored.png"), dpi=300)


#violin of scanpy_score and scanpy_score normalized
adata_auc.obs['THRB_score'] = adata_auc[:,"THRB(+)"].X.toarray().flatten()
adata_auc.obs['comparison'] = adata_auc_thrb.obs['comparison']
print(adata_auc.X[0:5,0:5])

adata = sc.read(os.path.join(path_data, "clustered_norm.h5ad"))
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

pairs = [
    ['promet_NT', 'nonpro_NT'],
    ['promet_AC', 'nonpro_AC'],
    ['promet_AC', 'promet_NT'],
    ['nonpro_AC', 'nonpro_NT'],
    ['nonpro_AC','promet_NT'],
    ['nonpro_NT','promet_AC']
]
order= ['nonpro_NT', 'nonpro_AC', 'promet_NT','promet_AC']
violin(
    df=adata_sc.obs,
    x='comparison',
    y='THRB_score',
    ax=ax,
    add_stats=True,
    pairs=pairs,
    linewidth=0.5
)

# Format the axis
plu.format_ax(ax, title='THRB scores', 
          xticks=adata_sc.obs['comparison'].cat.categories, ylabel='Score',reduced_spines=True)
ax.set_title('THRB scores', fontsize=16, fontweight='bold')
ax.set_xlabel('Condition', fontsize=14)
ax.set_ylabel('Score', fontsize=14)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, 'THRB_violin_scanpyscores_normalized_all_tests.png'), dpi=400)


#volcano plot
adata_auc.write(os.path.join(path_data,"clustered_scanpy_no_norm.h5ad"))

#create column percent_cells
regulon_expr_matrix = adata_sc.X > 0
percent_cells = np.sum(regulon_expr_matrix, axis=0) / regulon_expr_matrix.shape[0] * 100
adata_sc.var['percent_cells'] = np.ravel(percent_cells)


#add highly_variable_features
highly_variable_features = adata.var['highly_variable_features']
regulon_to_gene_mapping = {}  
for regulon_name in adata_sc.var.index:
    gene_name = regulon_name  
    regulon_to_gene_mapping[regulon_name] = gene_name

adata_sc.var['highly_variable_features'] = [
    highly_variable_features.get(regulon_to_gene_mapping[regulon], False)  
    for regulon in adata_sc.var.index
]

# add layer RAW 
#adata_tmp.layers['raw'] = adata_tmp.X

#add 'mean' and 'var' columns
T = adata_sc.X

if sparse.issparse(T):
    T = T.toarray()

adata_sc.X = T
mean_values = np.array(adata_sc.X.mean(axis=0))
variance_values = np.array(adata_sc.X.var(axis=0))

adata_sc.var['mean'] = mean_values
adata_sc.var['var'] = variance_values


# Check if the matrix is dense, if so, convert it to sparse CSR matrix
if isinstance(T, np.ndarray):
    print("Converting dense matrix to sparse CSR format...")
    adata_sc.X = csr_matrix(T)

if adata_sc.X.dtype != 'float32':
    adata_sc.X = adata_sc.X.astype('float32')


## DE ##

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata_sc, path_data, contrasts_name='paep_contrasts')

# Here we go
D = Dist_features(adata_sc, contrasts, jobs=jobs)
D.select_genes()
for k in D.jobs:
    for x in D.jobs[k]:
        if x['model'] == 'wilcoxon':
            job_key = '|'.join([k, x['features'], x['model']])
            de_results, gene_set_dict = D.compute_DE(contrast_key=k, which='perc_1_no_miribo')
            D.Results.add_job_results(de_results, gene_set_dict, job_key=job_key)

dfs = []
categories= ['nonpro_AC_vs_NT', 'promet_AC_vs_NT', 'promet_AC', 'pro_nonpro_AC', 'promet_NT']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, "Degs_regulon_scanpy_score.csv"))

promet=all_degs[all_degs['comparison']== 'promet_AC_vs_NT0_vs_promet_AC_vs_NT1']
promet['effect_size']
promet.to_csv(os.path.join(path_results,'DEgs_promet_AC_vs_promet_NT.csv'))
all_degs['comparison'].unique()



# #def compute_GSEA_for_single_contrast(all_degs, contrast, covariate='effect_size',
#                                      collection='GO_Biological_Process_2023', n_out=50):
#     """
#     Perform Gene Set Enrichment Analysis (GSEA) for a single contrast using effect_size (logFC).
    
#     Parameters:
#         - regulon_results: DataFrame containing regulon results with multiple contrasts
#         - contrast: Contrast to process (e.g., 'nonpro_AC_vs_NT')
#         - covariate: Metric to rank regulons by (e.g., 'effect_size')
#         - by: Sorting criteria for GSEA results (e.g., 'Adjusted P-value')
#         - collection: Gene set collection for GSEA (e.g., 'GO_Biological_Process_2021')
#         - n_out: Number of top results to return
    
#     Returns:
#         - GSEA results for the contrast
#     """
all_degs['comparison'].unique()
n_out=50
collection='GO_Biological_Process_2023'    #MSigDB_Hallmark_2020
contrast='promet_AC0_vs_promet_AC1'
contrast_results = all_degs[all_degs['comparison'] == contrast]
covariate='effect_size'
#if covariate not in contrast_results.columns:
    #raise ValueError(f"Covariate '{covariate}' not found in the regulon results for contrast '{contrast}'.")

ranked_regulon_list = contrast_results[['regulon', covariate]].dropna()
ranked_regulon_list = ranked_regulon_list.sort_values(covariate, ascending=False)  
ranked_regulon_list = ranked_regulon_list.set_index('regulon')[covariate]

#Perform GSEA
results = prerank(
    rnk=ranked_regulon_list,
    gene_sets=[collection],
    threads=cpu_count(),
    min_size=15,
    max_size=500,
    permutation_num=200,
    outdir=None,
    seed=1234,
    verbose=True,
)

df = results.res2d.loc[:, ['Term', 'ES', 'NES', 'FDR q-val','Lead_genes']].rename(columns={'FDR q-val': 'Adjusted P-value'})
df['Term'] = df['Term'].map(lambda x: x.split('__')[1] if '__' in x else x)
df = df.set_index('Term')
df.to_excel(os.path.join(path_results,'Gsea_promet_AC_vs_nonpromet_AC_GO.xlsx'))

#Viz GSEA 
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

plu.stem_plot(
    df[df['Adjusted P-value'] < 0.1][['ES', 'NES', 'Adjusted P-value']].sort_values('NES', ascending=False).head(25),
    'NES',
    ax=ax
)
plu.format_ax(ax, xlabel='NES')

fig.suptitle(f'GSEA: Prometastatic AC vs non-Prometastatic AC')
fig.tight_layout()
fig.savefig(os.path.join(path_results,"gsea_Prometastatic_AC_vs_nonprometastatic_AC_GO_filtered.png"),dpi=300)



#volcano plot
def get_genes_to_annotate_plot(df_, evidence, effect_size, n):
    df_ = df_[(df_['type'] != 'other')].copy()  
    df_['score'] = np.abs(df_[effect_size]) * df_[evidence]  
    return df_.sort_values('score', ascending=False).head(n).index.tolist()


def volcano_plot_plot(
    df, effect_size='effect_size', evidence='evidence',
    t_logFC=1, t_FDR=.1, n=10, title=None, xlim=(-8,8), max_distance=0.5, pseudocount=0,
    figsize=(5,5), annotate=False, s=5, lines=False
    ):
    """
    Volcano plot
    """    

    df_ = df.copy()    
    choices = [
        (df_[effect_size] >= t_logFC) & (df_[evidence] <= t_FDR),
        (df_[effect_size] <= -t_logFC) & (df_[evidence] <= t_FDR),
    ]
    df_['type'] = np.select(choices, ['up', 'down'], default='other')
    df_['to_annotate'] = False
    genes_to_annotate = get_genes_to_annotate_plot(df_, evidence, effect_size, n)
    df_.loc[genes_to_annotate, 'to_annotate'] = True
    df_[evidence] = -np.log10(df_[evidence]+pseudocount)

    fig, ax = plt.subplots(figsize=figsize)
    plu.scatter(df_.query('type == "other"'), effect_size, evidence,  c='darkgrey', s=s, ax=ax)
    plu.scatter(df_.query('type == "up"'), effect_size, evidence,  c='red', s=s*2, ax=ax)
    plu.scatter(df_.query('type == "down"'), effect_size, evidence,  c='b', s=s*2, ax=ax)

    ax.set(xlim=xlim)

    if lines:
        ax.vlines(1, df_[evidence].min(), df_[evidence].max(), colors='r')
        ax.vlines(-1, df_[evidence].min(), df_[evidence].max(), colors='b')
        ax.hlines(-np.log10(0.1), xlim[0], xlim[1], colors='k')

    plu.format_ax(ax, title=title, xlabel=f'log2FC', ylabel=f'-log10(FDR)')
    ax.spines[['top', 'right']].set_visible(False)

    if annotate:
        plu.ta.allocate_text(
            fig, ax, 
            df_.loc[lambda x: x['to_annotate']][effect_size],
            df_.loc[lambda x: x['to_annotate']][evidence],
            df_.loc[lambda x: x['to_annotate']].index,
            x_scatter=df_[effect_size], y_scatter=df_[evidence], 
            linecolor='black', textsize=8, 
            max_distance=max_distance, linewidth=0.5, nbr_candidates=100
        )

    return fig


all_degs['comparison'].unique()
prometnt=all_degs[all_degs['label']=='promet_AC']
all_degs= pd.read_csv(os.path.join(path_data, 'Degs_regulon_scanpy_score.csv'), index_col=0)
df_fixed= all_degs.copy()
df_fixed['evidence'] = all_degs['evidence'].replace(0, 1e-50)
df_fixed['evidence'] = df_fixed['evidence'].clip(lower=1e-50)
fig = volcano_plot_plot(
    df_fixed.query('comparison=="promet_AC_vs_NT0_vs_promet_AC_vs_NT1"'), effect_size='effect_size', evidence='evidence',
    n=30, annotate=True, xlim=(-2.5, 2.5),pseudocount=1e-50, title = "Prometastatic NT vs non-prometastatic NT"
)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results,"volcano_scanpy_score_norm_prometNT_vs_nonNT.png"),dpi=300)


#thrb genes
thrb_genes = d_reg['THRB(+)']['gene_set']

#umap single genes
adata_mp=adata.copy()
adata_mp.obs['comparison']= adata_auc_thrb.obs['comparison']
for gene in thrb_genes:
    plt.figure(figsize=(6, 6))  
    sc.pl.umap(adata_mp, color=gene, cmap='viridis', vmin=0,show=False)
    filename=f"umap_{gene}.png"
    plt.savefig(os.path.join(path_results,filename), dpi=300, bbox_inches='tight')

#mean expression genes belonging to thrb regulon in adata 
gene_means= {}

for gene in thrb_genes:
    if gene in adata_mp.var_names:
        gene_expr= adata_mp[:,gene].X.toarray().flatten()

        temp_df=pd.DataFrame({
            'expression': gene_expr,
            'comparison': adata_mp.obs['comparison'].values
        })

        mean_values= temp_df.groupby('comparison')['expression'].mean()
        gene_means[gene] = mean_values

#violin plot of each gene of thrb regulon
#for gene in thrb_genes:
    gene=['ANKS1A']
    adata_sc.obs[f'{gene}_score'] = adata[:, gene].X.toarray().flatten() 
    order= ['nonpro_NT', 'nonpro_AC', 'promet_NT','promet_AC']
    adata_sc.obs['comparison'] = pd.Categorical(
        adata_sc.obs['comparison'],
        categories=order,
        ordered=True
    )
    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111)

    pairs = [
        ['promet_NT', 'nonpro_NT'],
        ['promet_AC', 'nonpro_AC'],
        ['promet_AC', 'promet_NT'],
        ['nonpro_AC', 'nonpro_NT'],
        ['nonpro_AC','promet_NT']
    ]
    violin(
        df=adata_sc.obs,
        x='comparison',
        y=f'{gene}_score',
        ax=ax,
        add_stats=True,
        pairs=pairs,
        x_order=order,
        linewidth=0.5
    )

    # Format the axis
    plu.format_ax(ax, title=f'{gene} expression', 
            xticks=adata_sc.obs['comparison'].cat.categories, ylabel='Score',reduced_spines=True)
    ax.set_title(f'{gene} expression', fontsize=18, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=16)
    ax.set_ylabel('Expression', fontsize=16)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    fig.tight_layout()
    fig.savefig(os.path.join(path_results,f'{gene}_violin.png'), dpi=400)



#Thrb regulatory region
motif_paep = reg_paep['THRB(+)']['motif_data']
gene_hmga1 = d_reg['HMGA1(+)']['gene_set']
origin_thrb= regulons[regulons['TF']== 'THRB']
origin_hmga1= regulons[regulons['TF']== 'HMGA1']
origin_hmga1.to_csv(os.path.join(path_results,"hmga1_origin.csv"), index=False)
origin_thrb.to_csv(os.path.join(path_results,"thrb_origin.csv"), index=False)
origin_thrb= pd.read_csv(os.path.join(path_results,"thrb_origin.csv"))

df_= score_df_500[score_df_500['motifs']== 'tfdimers__MD00459']
targets_values = regulons[['TF','TargetGenes']]
target_values_thrb= targets_values[targets_values['TF']== 'THRB']

for x in origin_hmga1['TargetGenes']:
    print(x)
motif_paep = reg_paep['THRB(+)']['motif_data']





#umap all regulons
regulons_to_plot = adata_sc.var_names.tolist()  
for reg in regulons_to_plot:
    adata_mp.obs[f'{reg}_score'] = adata_mp[:,reg].X.toarray().flatten()

#Viz
plt.figure(figsize=(6, 6))  
sc.pl.umap(
    adata_mp,
    color=regulons_to_plot[7:14], 
    cmap='viridis',
    ncols=3
)
plt.show()
plt.savefig(os.path.join(path_results,"umap_comparison_regulons_all.png"), dpi=300, bbox_inches='tight')

     
#heatmap of regulon activity in contrasts
adata_sc.X = adata_sc.X.toarray()
df_ =pd.DataFrame(adata_sc.X, columns= adata_sc.var_names, index=adata_sc.obs_names)
df_['comparison']= adata_auc_thrb.obs['comparison']
mean_act=df_.groupby('comparison').mean().T
mean_act.loc['THRB']
reg_diff = mean_act.max(axis=1) - mean_act.min(axis=1)
selected_reg= reg_diff[reg_diff >= 0.7]
mean_act=mean_act.loc[selected_reg.index]
order=['nonpro_NT','nonpro_AC','promet_NT','promet_AC']
mean_act=mean_act[order]
plt.figure(figsize=(10,12))
plu.plot_heatmap(mean_act, vmin=-1, vmax=0.75, x_names_size=16, y_names_size=14)
plt.tight_layout()
plt.savefig(os.path.join(path_results,"heatmap_regulons_activity.png"), dpi=300, bbox_inches='tight')


#Jaccard similarity regulons (gene overlap)
d_reg['gene_set']= d_reg['gene_set'].apply(lambda x: set(x.split(',')))
d_reg['regulon']= d_reg['regulon'].str.replace('(+)', '', regex=False)
d_reg= d_reg[d_reg['regulon'].isin(selected_reg.index)]

def jaccardindex(set1,set2):
    intersection_size=len(set1 & set2)
    union_size= len(set1|set2)
    return intersection_size/union_size if union_size !=0 else 0

JI_matrix={}

for i, row_i in d_reg.iterrows():
    for j , row_j in d_reg.iterrows():
        if i >= j:
            continue
        regulon_i = row_i['regulon']
        regulon_j = row_j['regulon']
        gene_set_i = row_i['gene_set']
        gene_set_j = row_j['gene_set']

        JI= jaccardindex(gene_set_i, gene_set_j)
        JI_matrix[(regulon_i , regulon_j)] = JI
        JI_matrix[(regulon_j , regulon_i)] = JI

reg= d_reg['regulon'].values
JI_df= pd.DataFrame(np.zeros((len(reg), len(reg))), columns=reg, index=reg)

for (reg_a, reg_b), sim in JI_matrix.items():
    JI_df.loc[reg_a, reg_b] = sim
    JI_df.loc[reg_b, reg_a] = sim 

np.fill_diagonal(JI_df.values, 1)

order_clustering= leaves_list(linkage(JI_df.values,method='average'))
ordered_JI = JI_df.values[np.ix_(order_clustering, order_clustering)]
ordered_regulons = JI_df.index[order_clustering]
ordered_JI_df = JI_df.loc[ordered_regulons, ordered_regulons]
ordered_JI_df = pd.DataFrame(ordered_JI, index=JI_df.index[order_clustering], columns=JI_df.columns[order_clustering])
vmin, vmax= 0, 0.3
fig, ax = plt.subplots(figsize=(13, 13))
plu.plot_heatmap(ordered_JI_df, palette='mako', ax=ax,
                 x_names=ordered_regulons, y_names=ordered_regulons, annot=True, 
                 annot_size=8, label='Jaccard Index', shrink=1, cb=False)
sns.heatmap(ordered_JI_df, ax=ax, xticklabels=ordered_regulons, yticklabels=ordered_regulons,
                        robust=True, cmap="mako", vmin=vmin, vmax=vmax, fmt='.2f',cbar_kws={'fraction':0.05, 'aspect':35, 'pad': 0.02})
fig.tight_layout()
plt.savefig(os.path.join(path_results,"Jaccard_regulons_activity.png"), dpi=500, bbox_inches='tight')


#rna-seq enrichment analysis (bulk)
bulk_expr_mat=pd.read_csv(os.path.join(path_data,'bulk_expr_matr.csv'),index_col=0)

sample_to_condition = {
    'PT_shSCR_1': 'PT_shSCR', 'PT_shSCR_2': 'PT_shSCR', 'PT_shSCR_3': 'PT_shSCR',
    'PT_shSCR_4': 'PT_shSCR', 'PT_shSCR_5': 'PT_shSCR',
    'PT_shPAEP1_1': 'PT_shPAEP', 'PT_shPAEP1_2': 'PT_shPAEP', 'PT_shPAEP1_3': 'PT_shPAEP',
    'PT_shPAEP1_4': 'PT_shPAEP', 'PT_shPAEP1_5': 'PT_shPAEP',
    'PT_shPAEP2_3': 'PT_shPAEP', 'PT_shPAEP2_4': 'PT_shPAEP'
}

d_paep=d_reg['THRB(+)']['gene_set']
gene_set={'regulon': d_paep}

#Run ssgsea (bulk)
results = ssgsea(data=bulk_expr_mat, gene_sets=gene_set,
                 outdir=None, permutation_num=0, no_plot=True)
scores = results.res2d.set_index('Name')
scores['condition'] = scores.index.to_series().map(sample_to_condition)

#Viz
fig, ax = plt.subplots(figsize=(4,4))
order = ['PT_shPAEP','PT_shSCR']
scores["NES"] = pd.to_numeric(scores["NES"], errors="coerce")
plu.box(
    scores, 
    x='condition', y='NES', ax=ax, color='grey', add_stats=True, 
    pairs=[['PT_shPAEP','PT_shSCR']], 
    x_order=order
)
plu.strip(scores, x='condition', y='NES', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Regulon activity in shPAEP vs shSCR', ylabel='NES', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_shPAEP_regulon_act.png'), dpi=300)


#mean gene_set expression x condition (bulk)
#no PAEP and each gene of the regulon
d_paep=d_paep.remove('PAEP')
sub_bulk_matr= bulk_expr_mat.loc[bulk_expr_mat.index.intersection(d_paep)].T
sub_bulk_matr.columns
sub_bulk_matr['mean']= sub_bulk_matr.mean(axis=1)
sub_bulk_matr['condition'] = sub_bulk_matr.index.to_series().map(sample_to_condition)

#Viz
fig, ax = plt.subplots(figsize=(5.3,5))
order = ['PT_shPAEP','PT_shSCR']

plu.box(
    sub_bulk_matr, 
    x='condition', y='mean', ax=ax, color='grey', add_stats=True, 
    pairs=[['PT_shPAEP','PT_shSCR']], 
    x_order=order
)
plu.strip(sub_bulk_matr, x='condition', y='mean', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of regulon genes in shPAEP vs shSCR', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_shPAEP_regulon_mean_expr.png'), dpi=300)

#single gene (vers)
d_paep.remove("RAB33A")
for gene in d_paep:
    sub_bulk_matr= bulk_expr_mat.loc[[gene]].T
    sub_bulk_matr.columns = ['expression']
    sub_bulk_matr['condition'] = sub_bulk_matr.index.to_series().map(sample_to_condition)

    #Viz
    fig, ax = plt.subplots(figsize=(4,4))
    order = ['PT_shPAEP','PT_shSCR']

    plu.box(
        sub_bulk_matr, 
        x='condition', y='expression', ax=ax, color='grey', add_stats=True, 
        pairs=[['PT_shPAEP','PT_shSCR']], 
        x_order=order
    )
    plu.strip(sub_bulk_matr, x='condition', y='expression', ax=ax, color='k', x_order=order)
    plu.format_ax(ax=ax, title=f'{gene}: Expression in shPAEP vs shSCR', ylabel='expression', rotx=90, reduced_spines=True)
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, f'Boxplot_shPAEP_regulon_mean_expr_{gene}.png'), dpi=300)



"""
Script to analyze pyscenic results
"""

import os
import base64
import zlib
import json
import umap
import anndata as ad
import pandas as pd
import scanpy as sc
import loompy as lp
from gseapy import prerank
import gc
from multiprocessing import cpu_count
from scipy import sparse
from scipy.sparse import csr_matrix
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *
from Cellula.dist_features._Gene_set import Gene_set
from Cellula.plotting._plotting import *
from Cellula.dist_features._signatures import scanpy_score, wot_zscore
from matplotlib.gridspec import GridSpec
from plotting_utils.plotting_base import *
from Cellula._utils import sanitize_neighbors
from BC_chemo_utils.utils import prep_df_for_dotplot
from BC_chemo_utils.tests import *
from BC_chemo_utils.plotting import *
from Cellula.clustering._clustering import *
from Cellula.preprocessing._pp import * 
from Cellula.preprocessing._neighbors import *
from Cellula.preprocessing._embeddings import *
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















## DE ##

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






#Reformat regulons 
regulons.columns = regulons.iloc[0]         
regulons = regulons.drop(index=[0])     
regulons.columns.values[0] = regulons.iloc[0, 0]  
regulons.columns.values[1] = regulons.iloc[0, 1]
regulons= regulons.drop(index=[1]).reset_index(drop=True)

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
df= pd.read_csv(os.path.join(path_results,'regulon_data.csv'))

#Contrasts tested regulons
d_contrast={}
reg_to_extract= ['THRB(+)','SMAD5(+)','NR2F1(+)']
for key,value in d_reg.items():
    if key in  reg_to_extract:
        d_contrast[key] = value

#N regulons with PAEP
reg_paep = {}
for tf, value in d_reg.items():
    gene_set = value.get('gene_set', [])
    if 'PAEP' in gene_set:
        reg_paep[tf] = {
            'gene_set': gene_set,
            'motif_data': value.get('motif_data', {})
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
    params = update_params(params, kwargs)

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
                _cmap = create_palette(df, by, palette=categorical_cmap)
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
        add_wilcox(df, x, y, pairs, ax, order=x_order)

    return ax

def add_wilcox(
    df: pd.DataFrame, 
    x: str, 
    y: str, 
    pairs: Iterable[Iterable[str]], 
    ax:  matplotlib.axes.Axes = None, 
    order: Iterable[str] = None
    ):
    """
    Add statistical annotations (basic tests from statannotations).
    """
    annotator = Annotator(ax, pairs, data=df, x=x, y=y, order=order)
    annotator.configure(
        test='Mann-Whitney', text_format='star', show_test_name=False,
        line_height=0.001, text_offset=3
    )
    annotator.apply_and_annotate()


def format_ax(
    ax: matplotlib.axes.Axes = None, 
    title: str = None, 
    xlabel: str = None, 
    ylabel: str = None, 
    xticks: Iterable[Any] = None, 
    yticks: Iterable[Any] = None, 
    rotx: float = 0, 
    roty: float = 0, 
    axis: bool = True,
    xlabel_size: float = None, 
    ylabel_size: float = None,
    xticks_size: float = None, 
    yticks_size: float = None,
    title_size: float = None, 
    log: bool = False, 
    reduced_spines: bool = False
    ) -> matplotlib.axes.Axes:
    """
    Format labels, ticks and stuff of an ax: matplotlib.axes.Axes object.
    """

    if log:
        ax.set_yscale('log')
    
    if title is not None:
        ax.set(title=title)
    
    if xlabel is not None:
        ax.set(xlabel=xlabel)
    
    if xlabel is not None:
        ax.set(ylabel=ylabel)

    if xticks is not None:
        ax.set_xticks([ i for i in range(len(xticks)) ])
        ax.set_xticklabels(xticks)
    if yticks is not None:
        ax.set_yticks([ i for i in range(len(yticks)) ])
        ax.set_yticklabels(yticks)

    if xticks_size is not None:
        ax.xaxis.set_tick_params(labelsize=xticks_size)
    if yticks_size is not None:
        ax.yaxis.set_tick_params(labelsize=yticks_size)

    if xlabel_size is not None:
        ax.xaxis.label.set_size(xlabel_size)
    if ylabel_size is not None:
        ax.yaxis.label.set_size(ylabel_size)

    ax.tick_params(axis='x', labelrotation = rotx)
    ax.tick_params(axis='y', labelrotation = roty)

    if title_size is not None:
        ax.set_title(title, fontdict={'fontsize': title_size})
    
    if reduced_spines:
        ax.spines[['right', 'top']].set_visible(False)
    
    if not axis:
        ax.axis('off')

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
format_ax(ax, title='THRB scores', 
          xticks=adata_auc_thrb_clean.obs['comparison'].cat.categories, ylabel='Score', reduced_spines=True)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, 'THRB_violin.png'), dpi=400)


adata_auc_thrb.var.drop('THRB', inplace=True)

#umap
adata_mp=adata.copy()
adata_mp.obs['THRB_score'] = adata_tmp[:,"THRB"].X.toarray().flatten()
adata_mp.obs['comparison']= adata_auc_thrb.obs['comparison']
regulon_to_plot=['THRB_score_promet_NT']
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
plt.savefig(os.path.join(path_results,"umap_THRB_promet_NT.png"), dpi=300, bbox_inches='tight')
plt.close()







#adata_tmp=adata_auc_thrb_clean.copy()
# adata_tmp= adata_auc_thrb_clean[adata_auc_thrb_clean.obs['THRB'] != 0]
# adata_tmp.obs
pairs=[
    ['promet_NT', 'nonpro_NT'],
    ['promet_AC', 'nonpro_AC'],
    ['promet_AC', 'promet_NT'],
    ['nonpro_AC', 'nonpro_NT'],
    ['nonpro_AC','promet_NT']
]
pair= [['nonpro_NT', 'promet_AC']]
#median 
for group1, group2 in pair:
    median1= adata_auc_thrb.obs.query("comparison == @group1")['THRB'].median()
    median2= adata_auc_thrb.obs.query("comparison == @group2")['THRB'].median()
    
print(f"{group1} median: {median1}| {group2} median: {median2}")




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
W = adata_tmp.X

if sparse.issparse(W):
    W = W.toarray()

X_zscored= (W - W.mean(axis=0))/ W.std(axis=0)

adata_tmp.X = csr_matrix(X_zscored)

X_dense= adata_tmp.X.toarray()

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
adata_tmp
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111)

pairs = [
    ['promet_NT', 'nonpro_NT'],
    ['promet_AC', 'nonpro_AC'],
    ['promet_AC', 'promet_NT'],
    ['nonpro_AC', 'nonpro_NT'],
    ['nonpro_AC','promet_NT']
]
order= ['nonpro_NT', 'nonpro_AC', 'promet_NT','promet_AC']
violin(
    df=adata.obs,
    x='comparison',
    y='THRB_score',
    ax=ax,
    add_stats=True,
    pairs=pairs,
    linewidth=0.5
)

# Format the axis
format_ax(ax, title='THRB scores', 
          xticks=adata.obs['comparison'].cat.categories, ylabel='Score',reduced_spines=True)
ax.set_title('THRB scores', fontsize=16, fontweight='bold')
ax.set_xlabel('Condition', fontsize=14)
ax.set_ylabel('Score', fontsize=14)
ax.tick_params(axis='x', labelsize=12)
ax.tick_params(axis='y', labelsize=12)
fig.tight_layout()
plt.show()
fig.savefig(os.path.join(path_results, 'THRB_violin_scanpyscores_normalized.png'), dpi=400)


#volcano plot
adata_auc.write(os.path.join(path_data,"clustered_scanpy_no_norm.h5ad"))
adata_tmp = sc.read(os.path.join(path_data,"clustered_scanpy_no_norm.h5ad"))
adata_tmp.var['regulon'] = adata_tmp.var_names
adata_tmp.var['regulon'] = adata_tmp.var_names.str.replace('(+)', '', regex=False)
adata_tmp.var_names = adata_tmp.var['regulon']
adata_tmp.obs['THRB_score'] = adata_tmp[:,"THRB"].X.toarray().flatten()
adata_tmp.obs['comparison'] = adata_auc_thrb.obs['comparison']

#create column percent_cells
regulon_expr_matrix = adata_tmp.X > 0
percent_cells = np.sum(regulon_expr_matrix, axis=0) / regulon_expr_matrix.shape[0] * 100
adata_tmp.var['percent_cells'] = np.ravel(percent_cells)


#add highly_variable_features
highly_variable_features = adata.var['highly_variable_features']
regulon_to_gene_mapping = {}  
for regulon_name in adata_tmp.var.index:
    gene_name = regulon_name  
    regulon_to_gene_mapping[regulon_name] = gene_name

adata_tmp.var['highly_variable_features'] = [
    highly_variable_features.get(regulon_to_gene_mapping[regulon], False)  
    for regulon in adata_tmp.var.index
]

# add layer RAW 
adata_tmp.layers['raw'] = adata_tmp.X

#add 'mean' and 'var' columns
T = adata_tmp.X

if sparse.issparse(T):
    T = T.toarray()

adata_tmp.X = T
mean_values = np.array(adata_tmp.X.mean(axis=0))
variance_values = np.array(adata_tmp.X.var(axis=0))

adata_tmp.var['mean'] = mean_values
adata_tmp.var['var'] = variance_values


# Check if the matrix is dense, if so, convert it to sparse CSR matrix
if isinstance(T, np.ndarray):
    print("Converting dense matrix to sparse CSR format...")
    adata_tmp.X = csr_matrix(T)

if adata_tmp.X.dtype != 'float32':
    adata_tmp.X = adata_tmp.X.astype('float32')


## DE ##

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata_tmp, path_data, contrasts_name='paep_contrasts')

# Here we go
D = Dist_features(adata_tmp, contrasts, jobs=jobs)
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
all_degs=pd.read_csv(os.path.join(path_data, "Degs_regulon_scanpy_score.csv"))

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
n_out=50
collection='MSigDB_Hallmark_2020'
contrast='promet_AC_vs_NT0_vs_promet_AC_vs_NT1'
contrast_results = all_degs[all_degs['comparison'] == contrast]
covariate='effect_size'
#if covariate not in contrast_results.columns:
    #raise ValueError(f"Covariate '{covariate}' not found in the regulon results for contrast '{contrast}'.")

ranked_regulon_list = contrast_results[['regulon', covariate]].dropna()
ranked_regulon_list = ranked_regulon_list.sort_values(covariate, ascending=False)  
ranked_regulon_list = ranked_regulon_list.set_index('regulon')[covariate]

# Perform GSEA using gseapy.prerank
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
df.to_csv(os.path.join(path_results,'Gsea_promet_AC_vs_nonpro_AC.csv'))

#Viz GSEA 
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

stem_plot(
    df[['ES', 'NES', 'Adjusted P-value']].sort_values('NES', ascending=False).head(25),
    'NES',
    ax=ax
)
format_ax(ax, title='GSEA', xlabel='NES')

fig.suptitle(f'GSEA: Prometastatic AC vs Prometastatic NT')
fig.tight_layout()
fig.savefig(os.path.join(path_results,"gsea_Prometastatic_AC_vs_prometastatic_NT_hallmark'.png"),dpi=300)








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
    scatter(df_.query('type == "other"'), effect_size, evidence,  c='darkgrey', s=s, ax=ax)
    scatter(df_.query('type == "up"'), effect_size, evidence,  c='red', s=s*2, ax=ax)
    scatter(df_.query('type == "down"'), effect_size, evidence,  c='b', s=s*2, ax=ax)

    ax.set(xlim=xlim)

    if lines:
        ax.vlines(1, df_[evidence].min(), df_[evidence].max(), colors='r')
        ax.vlines(-1, df_[evidence].min(), df_[evidence].max(), colors='b')
        ax.hlines(-np.log10(0.1), xlim[0], xlim[1], colors='k')

    format_ax(ax, title=title, xlabel=f'log2FC', ylabel=f'-log10(FDR)')
    ax.spines[['top', 'right']].set_visible(False)

    if annotate:
        ta.allocate_text(
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
thrb_genes = reg_paep['THRB(+)']['gene_set']

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
for gene in thrb_genes:
    adata_mp.obs[f'{gene}_score'] = adata_mp[:, gene].X.toarray().flatten() 
    order= ['nonpro_NT', 'nonpro_AC', 'promet_NT','promet_AC']
    adata_mp.obs['comparison'] = pd.Categorical(
        adata_mp.obs['comparison'],
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
        df=adata_mp.obs,
        x='comparison',
        y=f'{gene}_score',
        ax=ax,
        add_stats=True,
        pairs=pairs,
        x_order=order,
        linewidth=0.5
    )

    # Format the axis
    format_ax(ax, title=f'{gene} scores', 
            xticks=adata_mp.obs['comparison'].cat.categories, ylabel='Score',reduced_spines=True)
    ax.set_title(f'{gene} scores', fontsize=16, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=14)
    ax.set_ylabel('Score', fontsize=14)
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, f'{gene}_violin.png'), dpi=400)





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



#boxplot sh conditions with PCA activity score for THRB
df_activity= pd.read_csv(os.path.join(path_results,"thrb_genes_expression_bulk.csv"),index_col=0)

fig, ax = plt.subplots(figsize=(4,4))
order = ['shSCR', 'shPAEP']
box(
    df_activity, 
    x='Condition', y='ActivityScore', ax=ax, c='grey', with_stats=True, 
    pairs=[['shSCR', 'shPAEP']], 
    order=order
)
strip(df_activity, 
    x='Condition', y='ActivityScore', ax=ax, c='k', order=order)
format_ax(ax=ax, title='n clones', ylabel='n', rotx=90)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'thrb_activity_shPAEP.png'), dpi=300)




     

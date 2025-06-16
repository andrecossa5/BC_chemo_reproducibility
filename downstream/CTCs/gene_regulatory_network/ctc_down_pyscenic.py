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
#adata_sc.write_h5ad(os.path.join(path_data,'clustered_norm.h5ad'))
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

#list TF of interest
tf_list = ["MYC","HIF1A","SREBF1", "SREBF2", "SREBP1","SREBP2","MLXIPL", 
           "CHREBP","PPARA","PPARD","PPARG","PPARB","TFAM","ATF4","FOXO1","TP53",
           "PTEN","MTOR", "MTORC1","HNF4A","HIF1","MTOR1"]
present_tf = [tf for tf in tf_list if tf in adata_sc.var_names]


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
    df_['is_present_tf'] = df_.index.isin(present_tf)
    df_['is_significant_tf'] = (
    df_['is_present_tf'] & 
    (np.abs(df_[effect_size]) >= t_logFC) & 
    (df_['evidence'] <= t_FDR)
)
    df_['to_annotate'] = False
    genes_to_annotate = set(get_genes_to_annotate_plot(df_, evidence, effect_size, n))
    genes_to_annotate |= set(df_.loc[df_['is_significant_tf']].index)  
    df_['to_annotate'] = df_.index.isin(genes_to_annotate)
    df_[evidence] = -np.log10(df_[evidence]+pseudocount)

    fig, ax = plt.subplots(figsize=figsize)
    plu.scatter(df_.query('type == "other"'), effect_size, evidence,  color='darkgrey', size=s, ax=ax)
    plu.scatter(df_.query('type == "up"'), effect_size, evidence,  color='red', size=s*2, ax=ax)
    plu.scatter(df_.query('type == "down"'), effect_size, evidence,  color='b', size=s*2, ax=ax)

    ax.set(xlim=xlim)

    if lines:
        ax.vlines(1, df_[evidence].min(), df_[evidence].max(), colors='r')
        ax.vlines(-1, df_[evidence].min(), df_[evidence].max(), colors='b')
        ax.hlines(-np.log10(0.1), xlim[0], xlim[1], colors='k')

    plu.format_ax(ax, title=title, xlabel=f'log2FC', ylabel=f'-log10(FDR)')
    ax.spines[['top', 'right']].set_visible(False)

    if annotate:
        highlight_mask = df_['to_annotate'] & df_['is_significant_tf']
        # normal_mask = df_['to_annotate'] & (~df_['is_significant_tf'])

        # ta.allocate_text(
        #     fig, ax,
        #     df_.loc[normal_mask, effect_size],
        #     df_.loc[normal_mask, evidence],
        #     df_.loc[normal_mask].index,
        #     x_scatter=df_[effect_size], y_scatter=df_[evidence],
        #     linecolor='black', textsize=8,
        #     max_distance=max_distance, linewidth=0.5, nbr_candidates=5
        # )

        ta.allocate_text(
            fig, ax,
            df_.loc[highlight_mask, effect_size],
            df_.loc[highlight_mask, evidence],
            df_.loc[highlight_mask].index,
            x_scatter=df_[effect_size], y_scatter=df_[evidence],
            linecolor='red', textsize=10, textcolor='red',
            max_distance=max_distance, linewidth=0.5, nbr_candidates=100
        )

    return fig


all_degs['comparison'].unique()
df_fixed= all_degs.copy()
df_fixed['evidence'] = all_degs['evidence'].replace(0, 1e-50)
df_fixed['evidence'] = df_fixed['evidence'].clip(lower=1e-50)
fig = volcano_plot_plot(
    df_fixed.query('comparison=="PT_m_vs_lung"'), effect_size='effect_size', evidence='evidence',
    n=30, annotate=True, xlim=(-2.5, 2.5),pseudocount=1e-50, title = "PT vs lung"
)
fig.tight_layout()
fig.savefig(os.path.join(path_results,"volcano_scanpy_score_norm_PT_vs_lung_tfselected.png"),dpi=300)
plt.close()

signi = df_fixed.query('comparison=="PT_m_vs_lung"')
matched_row= signi.loc[signi.index.intersection(present_tf)]

#violin plot 
for tf in present_tf:
    adata_sc.obs[f'{tf}_score'] = adata_sc[:,tf].X.toarray().flatten()

#violin function
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


plu.set_rcParams()

# Optional override for larger fonts
plt.rcParams.update({
    'font.size': 14,
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12
})
#Here we go
order= ['PT','lung','CTC']
for tf in present_tf:
    for col in adata_sc.obs.columns:
        if col == f'{tf}_score':
            fig = plt.figure(figsize=(10, 6))
            ax = fig.add_subplot(111)

            pairs = [
                ['PT', 'CTC'],
                ['lung', 'CTC'],
                ['PT', 'lung']
            ]

            violin(
                df=adata_sc.obs,
                x='origin',
                y=f'{tf}_score',   
                ax=ax,
                add_stats=True,
                pairs=pairs,
                x_order=order,
                linewidth=0.5
            )


            plu.format_ax(ax, title=f'{tf} scores', 
            xticks=order, ylabel='Score', reduced_spines=True)
            fig.tight_layout()
            fig.savefig(os.path.join(path_results, f'{tf}_violin.png'), dpi=400)


#UMAPS
for tf in present_tf:
    adata.obs[f'{tf}_score'] = adata_sc[:,tf].X.toarray().flatten()

#here we go
fig, ax = plt.subplots(figsize=(5, 5.5))
sc.pl.umap(adata, color='origin', cmap='viridis', vmin=0, show=False,size= 5, ax=ax)
plt.savefig(os.path.join(path_results,"umap_origin.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

for tf in present_tf:
    for col in adata.obs.columns:
        if col ==f'{tf}_score':
            fig, ax= plt.subplots(figsize=(5, 5.5))  
            sc.pl.umap(adata, color=f'{tf}_score', cmap='viridis', vmin=-3, vmax=2,show=False,size=5,ax=ax)
            plt.show()
            plt.savefig(os.path.join(path_results,f"umap_{tf}.png"), dpi=300, bbox_inches='tight')
            plt.close()



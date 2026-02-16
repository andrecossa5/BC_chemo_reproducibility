"""
Script to analyze pyscenic results CTC
Scoring used scanpy_score normalized
"""

import os
import base64
import zlib
import json
import umap
import numpy as np
import matplotlib.axes
from gseapy import ssgsea, prerank, enrichr
from multiprocessing import cpu_count
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
from BC_chemo_utils.tests import fastGSEA
from BC_chemo_utils.utils import prep_df_for_dotplot
from BC_chemo_utils.tests import *
matplotlib.use('macOSX')


#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility"
path_data= os.path.join(path_main,"data", "CTCs","grn","resources")
path_results= os.path.join(path_main, "results", "CTCs", "gene_reg_net","paper_ctc")

#Data
df_motif= pd.read_csv(os.path.join(path_main, "data","CTCs", "grn","resources", "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"),sep="\t")
regulons= pd.read_csv(os.path.join(path_data, "ctc_regulons.csv"))
adata= sc.read_h5ad(os.path.join(path_data, "clustered.h5ad"))
adata_sc = sc.read_h5ad(os.path.join(path_data, "clustered_norm.h5ad"))
lf= lp.connect(os.path.join(path_data, "ctc_aucell.loom"), mode="r+", validate=False)
auc_mtx= pd.DataFrame(lf.ca.RegulonsAUC, index=lf.ca.CellID)
metabolic_regulons= pd.read_csv(os.path.join(path_results, 'metabolic_regulons_filtered_msigdb.csv'))
# rank_df_10 = pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
# rank_df_500 = pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
# score_df_10=pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"))
# score_df_500= pd.read_feather(os.path.join(path_main, "data","CTCs", "grn","resources", "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"))

#filter auc_mtx for metabolism-related regulons
auc_mtx = auc_mtx[[col for col in auc_mtx.columns if col in metabolic_regulons['regulon'].values]]
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

#put all regulons as highly_variable_features
adata_auc.var['highly_variable_features'] = adata_auc.var['regulon']

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
plt.savefig(os.path.join(path_results, "distribution_aucell_metabolic.png"),dpi=300)


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
df = df[df['regulon'].isin(metabolic_regulons['regulon'].values)]
d_reg = {tf: d_reg[tf] for tf in metabolic_regulons['regulon'].values if tf in d_reg}
#save d_reg file

df.to_csv(os.path.join(path_results,'ctc_regulon.csv'), index=False)

#re-score values 
regulon_name = list(d_reg.keys())
len(regulon_name)
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
plt.savefig(os.path.join(path_results, "distribution_scanpyscores_metabolic.png"), dpi=300)

#normalized scanpy score 
Z = adata_sc.X
if sparse.issparse(Z):
    Z = Z.toarray()

X_zscored= (Z - Z.mean(axis=0))/ Z.std(axis=0)

adata_sc.X = csr_matrix(X_zscored)
adata_sc.write_h5ad(os.path.join(path_data,'clustered_norm_met_msig.h5ad'))
adata_sc= sc.read_h5ad(os.path.join(path_data,'clustered_norm_met_msig.h5ad'))
X_dense= adata_sc.X.toarray()

plt.figure(figsize=(8,5))
plt.hist(X_dense.flatten(), bins=50, color='steelblue', edgecolor= 'k')
plt.xlabel("Score")
plt.ylabel("Frequency")
plt.title("scanpy distribution (Z-score normalized)")
plt.tight_layout()
plt.savefig(os.path.join(path_results, "distribution_scanpy_norm_metabolic.png"), dpi=300)

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
categories= ['pt_ctc','lung_ctc','pt_lung','ctc_pt/lung']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, 'PT_CTC_DEGs__nozscore_metabolic_msig.csv'))

#list TF of interest
# tf_list = ["MYC","HIF1A","SREBF1", "SREBF2", "SREBP1","SREBP2","MLXIPL", 
#            "CHREBP","PPARA","PPARD","PPARG","PPARB","TFAM","ATF4","FOXO1","TP53",
#            "PTEN","MTOR", "MTORC1","HNF4A","HIF1","MTOR1"]
# present_tf = [tf for tf in tf_list if tf in adata_sc.var_names]


# #volcano plot 
def get_genes_to_annotate_plot(df_, evidence, effect_size, n):
    df_ = df_[(df_['type'] != 'other')].copy()  
    df_['score'] = np.abs(df_[effect_size]) * df_[evidence]  
    return df_.sort_values('score', ascending=False).head(n).index.tolist()


def volcano_plot_plot_classic(df, effect_size='effect_size', evidence='evidence', t_logFC=1, t_FDR=0.001, n=10, title=None, xlim=(-8,8), figsize=(5,5), annotate=False):
    df_ = df.copy()
    df_['type'] = np.select([(df_[effect_size]>=t_logFC)&(df_[evidence]<=t_FDR), (df_[effect_size]<=-t_logFC)&(df_[evidence]<=t_FDR)], ['up','down'], 'other')
    df_['annot'] = df_.index.isin(get_genes_to_annotate_plot(df_.query('type=="up"'), evidence, effect_size, n))
    df_[evidence] = -np.log10(df_[evidence]+1e-50)
    fig, ax = plt.subplots(figsize=figsize)
    for dt, c in [('other','grey'),('up','red'),('down','blue')]: plu.scatter(df_.query(f'type=="{dt}"'), effect_size, evidence, color=c, size=5 if dt=='other' else 10, ax=ax)
    ax.set(xlim=xlim); plu.format_ax(ax, title=title, xlabel='log2FC', ylabel='-log10(FDR)')
    if annotate: ta.allocate_text(fig, ax, df_.loc[df_['annot'], effect_size], df_.loc[df_['annot'], evidence], df_.loc[df_['annot']].index, x_scatter=df_[effect_size], y_scatter=df_[evidence])
    return fig

all_degs= pd.read_csv(os.path.join(path_data, 'PT_CTC_DEGs_metabolic_msig.csv')).set_index('regulon')
df_fixed= all_degs.copy()
df_fixed['evidence'] = all_degs['evidence'].replace(0, 1e-50)
df_fixed['evidence'] = df_fixed['evidence'].clip(lower=1e-50)

def get_comparison_title(label, direction):
    """Convert label and direction to biological comparison title"""
    label_map = {
        'pt_ctc': {'g0_vs_g1': 'PT vs CTC', 'g1_vs_g0': 'CTC vs PT'},
        'lung_ctc': {'g0_vs_g1': 'lung vs CTC', 'g1_vs_g0': 'CTC vs lung'}, 
        'pt_lung': {'g0_vs_g1': 'PT vs lung', 'g1_vs_g0': 'lung vs PT'},
        'ctc_pt/lung': {'g0_vs_g1': 'CTC vs PT/lung', 'g1_vs_g0': 'PT/lung vs CTC'}
    }
    return label_map.get(label, {}).get(direction, f"{label} ({direction})")

for label in df_fixed['label'].unique():
    df_label = df_fixed.query(f'label=="{label}"')
    for direction in ['g0_vs_g1', 'g1_vs_g0']:
        df_temp = df_label.query(f'comparison=="{direction}"')
        if len(df_temp) > 0:
            title = get_comparison_title(label, direction)
            filename = f"volcano_{label.replace('/', '_')}_{direction}_msig_prova.png"
            fig = volcano_plot_plot_classic(df_temp, effect_size='effect_size', evidence='evidence', n=30, annotate=True, xlim=(-2.5, 2.5), title=title)
            fig.tight_layout()
            fig.savefig(os.path.join(path_results, filename), dpi=300)
            plt.close()
            print(f"Saved {filename}")

#ORA for metabolic regulon classification 
def ORA_regulon_classification(d_reg, collection='MSigDB_Hallmark_2020', organism='human', n_terms=50):
    """ORA for regulon classification. Returns top N enriched terms per regulon (adjusted p-value <= 0.1)."""
    results = []
    for tf, genes in [(t, d['gene_set']) for t, d in d_reg.items() if len(d['gene_set']) >= 2]:
        try:
            df = enrichr(gene_list=genes, gene_sets=[collection], organism=organism, outdir=None).results
            df_filtered = df[df['Adjusted P-value'] <= 0.1]
            terms = ' | '.join(df_filtered['Term'].head(n_terms).tolist()) if len(df_filtered) > 0 else 'None'
            results.append({'regulon': tf, 'n_genes': len(genes), 'terms': terms})
        except: pass
    return pd.DataFrame(results)

# GO-biological ORA
ora_classification_go = ORA_regulon_classification(d_reg, collection='GO_Biological_Process_2025', n_terms=50)
ora_classification_go.to_csv(os.path.join(path_results, 'regulon_classification_GO.csv'), index=False)

# KEGG ORA
ora_classification_kegg = ORA_regulon_classification(d_reg, collection='KEGG_2021_Human', n_terms=50)
ora_classification_kegg.to_csv(os.path.join(path_results, 'regulon_classification_KEGG.csv'), index=False)

ora_classification_msig = ORA_regulon_classification(d_reg, collection='MSigDB_Hallmark_2020', n_terms=50)
ora_classification_msig.to_csv(os.path.join(path_results, 'regulon_classification_MSIGDB.csv'), index=False)

#filter_ora_classification_msigdb
count_terms = lambda ts, p: sum(1 for t in (ts.split('|') if isinstance(ts, str) and ts != 'None' else []) if re.search(p, t, re.IGNORECASE))
ora_classification_msig = pd.read_csv(os.path.join(path_results, 'regulon_classification_MSIGDB.csv'))
msig_kw = r"(?:metabol|glycol|fatty\s*acid|fattyacid|CoA|TCA|citric|krebs|catabo|anabol|lipid|" \
        r"amino\s*acid|nucleotide|purine|pyrimidine|gluconeogenesis|pentose\s*phosphate|" \
        r"oxidative\s*phosphorylation|electron\s*transport|beta-oxidation|NADH|FADH2|biosynth|hypoxia|HIF)"
ora_classification_msig['n_met_msigdb'] = ora_classification_msig['terms'].apply(lambda x: count_terms(x, msig_kw))
metabolic_regulons = ora_classification_msig[ora_classification_msig['n_met_msigdb'] > 0]
metabolic_regulons.to_csv(os.path.join(path_results, 'metabolic_regulons_filtered_msigdb.csv'), index=False)

# Plot top 15 regulons enriched by origin
X = adata_sc.X.toarray() if sparse.issparse(adata_sc.X) else adata_sc.X
df_activity = pd.DataFrame(X, columns=adata_sc.var_names, index=adata_sc.obs_names)
df_activity['origin'] = adata_sc.obs['origin'].values
mean_activity = df_activity.groupby('origin')[adata_sc.var_names].median()
fig, axs = plt.subplots(1, 3, figsize=(16, 6))
for ax, origin in zip(axs, ['PT', 'lung', 'CTC']):
    top_regs = mean_activity.loc[origin].nlargest(15)
    top_regs.plot(kind='barh', ax=ax, color='steelblue')
    plu.format_ax(ax=ax, title=f'Top 15 regulons in {origin}', xlabel='Mean activity', reduced_spines=True)
    ax.invert_yaxis()

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'origin_enrichment.png'), dpi=300)

"""
Script to analyze senescence signatures
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
path_data= os.path.join(path_main,"data", "MDA")
#path_adata = os.path.join("/Users/ieo7295/Desktop/tests/cell/senescence_paep/data/SCT")
path_results= os.path.join(path_main, "results", "MDA")

#Data
adata= sc.read_h5ad(os.path.join(path_data, "clustered.h5ad"))
adata_auc_thrb= sc.read(os.path.join(path_data,"clustered_thrb.h5ad"))
sign_1 = pd.read_csv(os.path.join(path_data, "REACTOME_CELLULAR_SENESCENCE.v2024.1.Hs.tsv"), sep='\t')
sign_2 = pd.read_csv(os.path.join(path_data, "SAUL_SEN_MAYO.v2024.1.Hs.tsv"), sep='\t')
sign_3 = pd.read_csv(os.path.join("/Users/ieo7295/Downloads/FRIDMAN_SENESCENCE_UP.v2024.1.Hs.tsv"),sep='\t')
sign_4 = ['IFIT1','IFIT3','OAS1','SYTL2','IFIT2','RTP4','OASL','N4BP2L1','IFITM1','AMY2B','PARP9','IFI44','XAF1','ISG15','CLIC5','IFI6','IFI27','RSAD2',
'NR4A2','BLNK','OAS2','CMPK2']

with open(os.path.join(path_data,'signatures.pickle'), 'rb') as f:
    data = pickle.load(f)

reactome = sign_1['REACTOME_CELLULAR_SENESCENCE'][16].split(',')
mayo = sign_2['SAUL_SEN_MAYO'][16].split(',')
friedman = sign_3['FRIDMAN_SENESCENCE_UP'][16].split(',')

reactome_present = []
for gene in reactome:
    if gene in adata.var_names:
        reactome_present.append(gene)
len(reactome_present)
mayo_present = []
for gene in mayo:
    if gene in adata.var_names:
        mayo_present.append(gene)
len(mayo_present)
friedman_present = []
for gene in friedman:
    if gene in adata.var_names:
        friedman_present.append(gene)
len(friedman_present)

TIS_present = []
for gene in sign_4:
    if gene in adata.var_names:
        TIS_present.append(gene)
len(TIS_present)

common_fridman_mayo = [x for x in friedman_present if x in mayo_present]
common_fridman_reactome =[x for x in friedman_present if x in reactome_present]
common_mayo_reactome = [x for x in mayo_present if x in reactome_present]

with open('/Users/ieo7295/Desktop/tests/cell/senescence_paep/reactome_present_prova.txt', 'w') as f:
    for gene in reactome_present:
        f.write(f"{gene.strip()};")


with open('/Users/ieo7295/Desktop/tests/cell/senescence_paep/mayo_present.txt', 'w') as f:
    for gene in mayo_present:
        f.write(f"mayo_present\t{gene.strip()}\n")

with open('/Users/ieo7295/Desktop/tests/cell/senescence_paep/friedman_present.txt', 'w') as f:
    for gene in friedman_present:
        f.write(f"fridman_present\t{gene.strip()}\n")

with open('/Users/ieo7295/Desktop/BC_chemo_reproducibility/data/MDA/TIS_present.txt', 'w') as f:
    for gene in TIS_present:
        f.write(f"TIS_present\t{gene.strip()}\n")

#add_column comparison
adata.obs['comparison'] = adata_auc_thrb.obs['comparison']

#prepare data for AC vs NT contrast
import numpy as np
adata_deg = adata.copy()
adata_deg.obs["treatment"] = np.where(
    adata_deg.obs["comparison"].str.contains("AC"), "AC", "NT"
)
adata_deg.write(os.path.join(path_data, "clustered_treatment.h5ad"))

#compute paep , mayo and hUSI median for table
paep_median = adata.obs.groupby("comparison")["PAEP"].median().sort_values(ascending=False)
adata.obs['mayo_present']= data['scores']['mayo_present']
mayo_median = adata.obs.groupby("comparison")["mayo_present"].median().sort_values(ascending=False)
husi_median = adata.obs.groupby("comparison")["hUSI"].median().sort_values(ascending=False)
adata.obs['TIS']= data['scores']['TIS_present']
tis_median = adata.obs.groupby("comparison")["TIS"].median().sort_values(ascending=False)
#subset adata for only cells of interest
adata_comparison = adata[adata.obs['comparison'].isin(['promet_AC', 'nonpro_AC', 'promet_NT', 'nonpro_NT'])].copy()
adata_comparison.write_h5ad(os.path.join(path_data, 'clustered_comparison.h5ad'))
adata_comparison = sc.read_h5ad(os.path.join(path_data, 'clustered_comparison.h5ad'))
del adata_comparison.obs['PAEP']
del adata_comparison.obs['hUSI']
del adata_comparison.obs['sen_binary']
del adata_comparison.obs['sen_gmm']
#test sasp_specific signature extracted from MAyo 
# sasp_sign= ['CCL2','CCL7','CXCL1','CXCL2','CXCL3','CXCL8','CXCL12','IL1A','IL1B','IL6',
# 'IL10','IL13','IL15','IL18','MMP1','MMP3','MMP9','SERPINE1','VEGFC','ANGPTL4','IGFBP2',
# 'IGFBP3','IGFBP4','TNFRSF1A','TGFA','SPP1','TIMP2']

# sasp_all = ['IL6', 'IL1A', 'CXCL8', 'CXCL10', 'CXCL12', 'TNF', 'TNFAIP2', 'TNFAIP3',
#  'CXCL14', 'CCL5', 'CXCL1', 'CCL20', 'AREG', 'IGFBP1', 'IGFBP2', 'IGFBP3',
#  'IGFBP4', 'IGFBP5', 'IGFBP6', 'IGFBP7', 'GDF15', 'EGF', 'FGF2', 'VEGFA',
#  'HGF', 'MMP1', 'MMP10', 'MMP13', 'MMP3', 'MMP9', 'SERPINE1', 'SERPINB2',
#  'PLAU', 'THBS1', 'SPARC', 'TIMP2', 'CDKN1A', 'CDKN2D', 'TP53', 'STAT1',
#  'STAT3', 'FOS', 'JUN']

# sasp_unique = list(set(sasp_sign + sasp_all))

#scanpy scoring 
# d_sig = {
#     'mayo': mayo_present,
#    #'fridman' : friedman_present,
#     #'reactome' : reactome_present,
#     'TIS': TIS_present,
#    # 'sasp_sign' : sasp_sign,
#    # 'sasp_unique': sasp_unique

# }

# for key in d_sig:
#     gene_set = d_sig[key]
#     score_name = f"{key}_scanpy"
#     gene_check = [g for g in gene_set if g in adata_comparison.var_names]
#     scores = scanpy_score(adata_comparison,gene_check)
#     adata_comparison.obs[f'{key}_scanpy'] = scores.values

#cellula_scored signature.py
signatures = data['scores']
adata.obs['fridman_present'] = signatures['friedman_present']
adata.obs['mayo_present'] = signatures['mayo_present']
adata.obs['reactome_present'] = signatures['reactome_present']
adata.obs['PAEP'] = adata[:,'PAEP'].X.toarray().flatten()
#adata.write_h5ad(os.path.join(path_data, 'clustered_new.h5ad'))

#paep-senescence correlation
paep_expr = adata[:,'PAEP'].X.toarray().flatten()

df = pd.DataFrame({
    'PAEP': paep_expr,
    'TIS': adata.obs['TIS'],
    #'hUSI': adata.obs['hUSI'],
    #'mayo_scanpy' : adata_comparison.obs['mayo_scanpy'],
    # 'fridman_present' : adata.obs['fridman_present'],
    # 'reactome_present': adata.obs['reactome_present'], #    'sasp_sign_score': adata.obs['sasp_sign_score'],'sasp:unique_score': adata.obs['sasp_unique_score']
    'comparison': adata.obs['comparison']
})
colors = {
    'promet_NT': 'blue',
    'nonpro_NT': 'orange',
}

signatures = ['TIS'] #'sasp_sign_score','mayo_present','fridman_present','reactome_present','sasp:unique_score' ,'hUSI'
df = df[df['comparison'].isin(['promet_NT','nonpro_NT'])]
df['comparison'] = df['comparison'].cat.remove_unused_categories()


for signature in signatures:
    corr_res = {}
    for group in colors.keys():
        sub_df = df[df['comparison'] == group]
        r, p = spearmanr(sub_df['PAEP'], sub_df[signature])
        corr_res[f"{group}_{signature}"] = f"{group}, {signature}: r = {r:.2f}, p= {p:.1e}"


    sns.set_theme(style='whitegrid')
    g = sns.lmplot(
        data=df, 
        x='PAEP',
        y=signature,
        hue='comparison',
        palette=colors,
        height=5,
        aspect=1.2,
        scatter_kws={'s': 10, 'alpha': 0.5},
        line_kws={'linewidth': 2},
        legend=False
    )

    ax = g.ax
    y_max = df[signature].max()
    y_min = df[signature].min()
    offsets = [0, -0.05 * (y_max - y_min)]

    for i, group in enumerate(colors.keys()):
        key = f"{group}_{signature}"
        text = corr_res.get(key, '')
        ax.text(
            x=0.05 * (df['PAEP'].max() - df['PAEP'].min()) + df['PAEP'].min(),
            y=y_max + offsets[i],
            s=text,
            color=colors[group],
            fontsize=10,
            weight='bold'
        )

    plt.title(f'Correlation between PAEP and {signature.replace("_present", "").capitalize()} Senescence Signature')
    plt.xlabel('PAEP Expression (log-normalized)')
    plt.ylabel(f'{signature.replace("_present", "").capitalize()} Signature Score')
    plt.tight_layout()

    filename = f'corr_paep_{signature}_{colors.keys()}.png'
    plt.savefig(os.path.join(path_results, filename), dpi=300)
    plt.close()


#Co-localization PAEP - MAYO
df = adata.obs.copy()
df['PAEP'] = adata[:, 'PAEP'].X.toarray().flatten()
df['TIS'] = adata.obs['TIS']
df['UMAP1'] = adata.obsm['X_umap'][:, 0]
df['UMAP2'] = adata.obsm['X_umap'][:, 1]
del adata.obs['PAEP']
sc.pl.umap(
    adata,
    color=['TIS'],
    cmap='viridis',  
    vmax=[0.20, 6],
    title='TIS Score and PAEP Co-localization',
    show=True
)



fig, axs = plt.subplots(1, 2, figsize=(10, 5))

# Mayo score with vmax
sc.pl.umap(
    adata,
    color='TIS',
    cmap='viridis',
    vmax=0.7,
    ax=axs[0],
    show=False,
)
del adata.obs['PAEP']
# PAEP expression
sc.pl.umap(
    adata,
    color='PAEP',
    cmap='viridis',
    vmax=2,
    ax=axs[1],
    show=True,
    title='PAEP Expression'
)

plt.tight_layout()
plt.savefig(os.path.join(path_results, 'umap_tis_paep'), dpi=300)



promet_df= adata[adata.obs['comparison']== 'promet_AC']
promet_df.obs['PAEP'] = promet_df[:, 'PAEP'].X.toarray().flatten()
fridman = promet_df.obs['fridman_present']
mayo= promet_df.obs['mayo_present']
paep = promet_df.obs['PAEP']
reactome = promet_df.obs['reactome_present']

df= promet_df.obs.copy()


plt.figure(figsize=(6,5))
sns.scatterplot(x=reactome, y=paep, alpha=0.4)
plt.xlabel('reactome Score')
plt.ylabel('PAEP Expression')
plt.title('PAEP vs reactome')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(path_results, 'scatter-mayonon_paep'), dpi=300)


promet_ac_clones = [
    "CGGAAGTCCATCCCCTCG", "CTGATAAGGCGTCGAGTT", "TTGGGGGTAAACAGGGTC",
    "AACAATCCCACAGGAAGC", "GCGACAGCGTCTCTTTCC", "CATGTCCAGCTGTGGTGG"
]

promet_nt_clones = [
    "TACGGCGCGTCCTGTGGC","GCTGCATAGGCTGTAGCG","CACAATGCACTGCGAACT",
    "CCCAGCGACTCCGCACTG","GTCGCTGTCCTGCTCCCG","GCTATAACAAATTACATC",
    "CATACAGGTGACCGTGGC","CCCTTGCTTCCACTGTCC","CCTTTAGCAAGCCCGAAA",
    "CACGATCCGTCGGCCCAG","ATGTATCAGGTTCGCGAA","CGTCTCTAAGCGTGAGTC",
    "TTCGGCCAGGCTTTCCGT","CGCTTTCACCCGTCGAGT","CTTGGGTACGCGTATCGG",
    "AGATGATACACGTGGGCT"
]


subset_mask = (
    (adata.obs['condition'] == "PT, treated") &
    (adata.obs['GBC'].isin(promet_ac_clones))
)
adata_subset = adata[subset_mask]

for clone in adata_subset.obs['GBC'].unique():
    clone_data = adata_subset[adata_subset.obs['GBC'] == clone]
    mayo = clone_data.obs['mayo_present']
    paep = clone_data[:, 'PAEP'].X.toarray().flatten()

    if len(mayo) >= 30:  
        r, p = spearmanr(mayo, paep)
        print(f"GBC clone {clone}: Spearman r = {r:.2f}, p = {p:.3e}, n = {len(mayo)}")
    else:
        print(f"GBC clone {clone}: not enough cells (n={len(mayo)})")


sub_df = (
    (adata.obs['condition']== 'PT, treated') &
    (adata.obs['GBC']== 'GCGACAGCGTCTCTTTCC')
)
an_sub_df = adata[sub_df]
paep = an_sub_df[:,'PAEP'].X.toarray().flatten()
mayo= an_sub_df.obs['mayo_present']
fridman= an_sub_df.obs['fridman_present']
plt.figure(figsize=(6,5))
sns.scatterplot(x=mayo, y=paep, alpha=0.4)
plt.xlabel('reactome Score')
plt.ylabel('PAEP Expression')
plt.title('PAEP vs reactome')
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(path_results, 'scatter-mayonon_paep'), dpi=300)


#hUSI scoring
import sys
sys.path.append('/Users/ieo7295/Desktop/HUSI')
from sc.hUSI import cal_hUSI, SSE_hUSI, GMM_hUSI
import numpy as np
import matplotlib.colors as clr
import rpy2.robjects.packages as rpackages
import scanpy as sc
from rpy2.robjects.packages import importr
mclust = importr('mclust')

utils = rpackages.importr('utils')
utils.install_packages('mclust')
hUSI = cal_hUSI(adata)
adata.obs['hUSI'] = hUSI
color_self = clr.LinearSegmentedColormap.from_list('pink_grey', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
sc.pl.tsne(adata,color='hUSI',save='Python_demo_hUSI.png',size=30,cmap = color_self)

sc.pl.umap(adata_comparison, color='hUSI', size=10,cmap=color_self, show=True)
plt.savefig(os.path.join(path_results,'umap_hUSI'), dpi=300, bbox_inches='tight')

SenClass_binary = SSE_hUSI(hUSI)
adata.obs['sen_binary']= SenClass_binary

SenClass = GMM_hUSI(hUSI)
adata.obs['sen_gmm'] = SenClass
adata.write(os.path.join(path_data, "clustered_hUSI.h5ad"))
adata= sc.read_h5ad(os.path.join(path_data, "clustered_hUSI.h5ad"))
adata.obs['PAEP'] = adata[:,'PAEP'].X.toarray().flatten()

hUSI_score = adata.obs['hUSI']

adata_score = adata.copy()
adata_score = adata_score[adata_score.obs['hUSI']> 0.75]
sc.pl.umap(adata, color='hUSI,PAEP', size=10,cmap=color_self, show=True)

fig, axs = plt.subplots(1, 2, figsize=(10, 5))
mask = adata.obs['comparison'] == 'nonpro_AC'
adata.obs['hUSI_nonpro_AC'] = adata.obs['hUSI'].where(mask)
sc.pl.umap(adata, color='hUSI_nonpro_AC', ax=axs[0],size=10,cmap=color_self, show=False)
adata.obs['PAEP_nonpro_AC'] = adata.obs['PAEP'].where(mask)
sc.pl.umap(adata, color='PAEP_nonpro_AC',ax=axs[1], size=10,cmap=color_self, show=False)
plt.tight_layout()
plt.savefig(os.path.join(path_results, 'umap_hUSI_paep_grey_nonproAC'), dpi=300)


#binarization of hUSI and mayo 
adata.obs['sen_gmm'].unique()
adata.obs['mayo_present_bin'] = (adata.obs['mayo_present'] > 0.5).astype(int)

cluster_means = adata.obs.groupby("sen_gmm")["hUSI"].mean()
print(cluster_means.sort_values(ascending=False))
adata.obs["hUSI_binary_gmm"] = (adata.obs["sen_gmm"] == 4).astype(int)

sen_cluster = cluster_means.idxmax()
adata.obs["hUSI_binary_gmm"] = (adata.obs["sen_gmm"] == sen_cluster).astype(int)

summary = adata.obs.groupby("comparison")["hUSI_binary_gmm"].agg(["mean", "sum", "count"])
summary.rename(columns={"mean": "%_senescent", "sum": "n_senescent", "count": "n_total"}, inplace=True)
summary["%_senescent"] *= 100
print(summary.sort_values("%_senescent", ascending=False))

adata.obs['TIS'] = data['scores']['TIS_present']
tis_scores = adata.obs['TIS'].copy()

tis_scores.to_csv(os.path.join(path_data, "tis_scores.csv"), index=True)
tis_scores= pd.read_csv(os.path.join(path_data, "tis_scores.csv"), index_col=0)
adata.obs['mayo']= data['scores']['mayo_present']

mayo_scores = adata.obs["mayo"].copy()
sen_mayo = GMM_hUSI(mayo_scores)
adata.obs["mayo_gmm"] = sen_mayo
cluster_means = adata.obs.groupby("mayo_gmm")["mayo"].mean()
print(cluster_means.sort_values(ascending=False))
sen_cluster = cluster_means.idxmax()
adata.obs["mayo_binary_gmm"] = (adata.obs["mayo_gmm"] == sen_cluster).astype(int)

adata.obs["mayo_binary_gmm"] = (adata.obs["mayo_gmm"] == 4).astype(int)
summary = adata.obs.groupby("comparison")["mayo_binary_gmm"].agg(["mean", "sum", "count"])
summary.rename(columns={"mean": "%_senescent", "sum": "n_senescent", "count": "n_total"}, inplace=True)
summary["%_senescent"] *= 100
print(summary.sort_values("%_senescent", ascending=False))

tis_score = adata.obs['TIS'].copy()
sen_tis = GMM_hUSI(tis_score)
adata.obs["tis_gmm"] = sen_tis
cluster_means = adata.obs.groupby("tis_gmm")["TIS"].mean()
print(cluster_means.sort_values(ascending=False))
sen_cluster = cluster_means.idxmax()
adata.obs["tis_binary_gmm"] = (adata.obs["tis_gmm"] == sen_cluster).astype(int)

adata.obs["tis_binary_gmm"] = (adata.obs["tis_gmm"] == 4).astype(int)
summary = adata.obs.groupby("comparison")["tis_binary_gmm"].agg(["mean", "sum", "count"])
summary.rename(columns={"mean": "%_senescent", "sum": "n_senescent", "count": "n_total"}, inplace=True)
summary["%_senescent"] *= 100
print(summary.sort_values("%_senescent", ascending=False))



#DE promet_AC vs nonpro_AC in AC_AC_PTs_2
subset_cells = adata.obs["sample"] == "AC_AC_PTs_2"  
adata_half = adata[subset_cells].copy()
adata_half.write(os.path.join(path_data, "AC_AC_PTs_2.h5ad"))
adata_half = sc.read_h5ad(os.path.join(path_data, "AC_AC_PTs_2.h5ad"))

# Prepare jobs and contrasts for DE analysis
jobs, contrasts = prep_jobs_contrasts(adata_deg, path_data, contrasts_name='pro_AC')

# Run DE using Wilcoxon model
D = Dist_features(adata_deg, contrasts, jobs=jobs)
D.select_genes()
for k in D.jobs:
    for x in D.jobs[k]:
        if x['model'] == 'wilcoxon':
            job_key = '|'.join([k, x['features'], x['model']])
            de_results, gene_set_dict = D.compute_DE(contrast_key=k, which='perc_1_no_miribo')
            D.Results.add_job_results(de_results, gene_set_dict, job_key=job_key)

dfs = []
categories= ['promet_AC']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
df_mayo = all_degs.loc[all_degs.index.isin(mayo_present)].copy()
df_mayo.index.name = "gene"
df_mayo = df_mayo.reset_index()
df_mayo= df_mayo.sort_values("effect_size", ascending=False)



def get_genes_to_annotate_plot(df_, evidence, effect_size, n):
    df_ = df_[(df_['type'] != 'other')].copy()  
    df_['score'] = np.abs(df_[effect_size]) * df_[evidence]  
    return df_.sort_values('score', ascending=False).head(n).index.tolist()


def volcano_plot_plot(
    df, effect_size='effect_size', evidence='evidence',
    t_logFC=0.5, t_FDR=.1, n=10, title=None, xlim=(-8,8), max_distance=0.5, pseudocount=0,
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

#volcano_plot
df_fixed= df_mayo.copy()
df_fixed['evidence'] = df_mayo['evidence'].replace(0, 1e-50)
df_fixed['evidence'] = df_fixed['evidence'].clip(lower=1e-50)
fig = volcano_plot_plot(
    df_fixed.query('comparison=="g0_vs_g1"'), effect_size='effect_size', evidence='evidence',
    n=30, annotate=True, xlim=(-2.5, 2.5),pseudocount=1e-50, title = "Promet AC vs non-promet AC (Mayo)",
)
fig.tight_layout()
fig.savefig(os.path.join(path_results,"volcano_all_prometAC_vs_nonAC.png"),dpi=300)

#hUSI and PAEP quantification
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import scipy.stats as st

# ----- Data -----
df = adata.obs[["TIS","PAEP","comparison"]]
df = df.dropna(subset=["comparison","PAEP","TIS"])


# ------------------------
# 3. Boxplots with condition-specific tertiles + Cliff's δ & Δmedian
# ------------------------
def tertiles_by_rank(s):
    pct = s.rank(pct=True, method="average")
    return pd.cut(pct, bins=[0, 1/3, 2/3, 1], labels=["Low","Mid","High"], include_lowest=True)

df["PAEP_bin_cond"] = df.groupby("comparison", observed=False)["PAEP"].transform(tertiles_by_rank)

def cliffs_delta(x, y):
    x = np.asarray(x); y = np.asarray(y)
    return (np.sum(x[:,None] > y) - np.sum(x[:,None] < y)) / (len(x)*len(y))

rows = []
for comp, sub in df.dropna(subset=["PAEP_bin_cond"]).groupby("comparison"):
    #low  = sub.loc[sub.PAEP_bin_cond=="Low","hUSI"]
    high = sub.loc[sub.PAEP_bin_cond=="High","TIS"]
    medium = sub.loc[sub.PAEP_bin_cond=="Mid","TIS"]
    p = st.mannwhitneyu(medium, high, alternative="two-sided").pvalue
    delta = cliffs_delta(high, medium)
    dmed = high.median() - medium.median()
    rows.append({"comparison": comp, "p_raw": p, "delta": delta, "d_median": dmed,
                 "n_medium": len(medium), "n_high": len(high)})

box_stats = pd.DataFrame(rows)
box_stats["p_adj"] = multipletests(box_stats["p_raw"], method="fdr_bh")[1]
print("\nBoxplot stats (cond-specific tertiles):")
print(box_stats)

# Plot
comps = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']                               #sorted(df["comparison"].unique())
fig, axes = plt.subplots(1, len(comps), figsize=(5*len(comps), 5), sharey=True)
axes = np.ravel([axes]) 

for ax, comp in zip(axes, comps):
    sub = df[df.comparison==comp]
    data = [sub.loc[sub.PAEP_bin_cond==lab,"TIS"] for lab in ["Low","Mid","High"]]
    bp = ax.boxplot(data, positions=[0,1,2], widths=0.6, patch_artist=True, showfliers=False)
    ax.set_xticks([0,1,2]); ax.set_xticklabels(["Low","Mid","High"], fontsize=12)
    if ax is axes[0]: ax.set_ylabel("TIS", fontsize=12)
    ax.set_title(comp, fontsize=14)

    row = box_stats[box_stats.comparison==comp].iloc[0]
    txt = f"Mid vs High p={row.p_adj:.2e}\nδ={row.delta:.2f}"
    ax.text(0.5, 0.95, txt, transform=ax.transAxes, ha="center", va="top", fontsize=12)

plt.tight_layout()
plt.savefig(os.path.join(path_results, 'boxplots_tertiles_medium_TIS'), dpi=300)

#To compare effect sizes between two conditions, use a bootstrap approach to estimate the difference in Cliff's δ between the two groups. 

def bootstrap_cliffs_delta_diff(df, group_col, value_col, condition_col, cond1, cond2, n_iter=1000):
    deltas1, deltas2 = [], []
    for _ in range(n_iter):
        samp1 = df[df[condition_col] == cond1].sample(frac=1, replace=True)
        samp2 = df[df[condition_col] == cond2].sample(frac=1, replace=True)

        low1 = samp1[samp1[group_col] == "Low"][value_col]
        high1 = samp1[samp1[group_col] == "High"][value_col]
        delta1 = cliffs_delta(high1, low1)

        low2 = samp2[samp2[group_col] == "Low"][value_col]
        high2 = samp2[samp2[group_col] == "High"][value_col]
        delta2 = cliffs_delta(high2, low2)

        deltas1.append(delta1)
        deltas2.append(delta2)

    delta_diff = np.array(deltas1) - np.array(deltas2)
    ci_low, ci_high = np.percentile(delta_diff, [2.5, 97.5])
    p_val = np.mean(delta_diff <= 0) if np.mean(delta_diff) > 0 else np.mean(delta_diff >= 0)
    return ci_low, ci_high, p_val

# Example usage:
ci_low, ci_high, p_val = bootstrap_cliffs_delta_diff(
    df, group_col="PAEP_bin_cond", value_col="hUSI",
    condition_col="comparison", cond1="nonpro_AC", cond2="nonpro_NT"
)

print(f"Δδ (NT - NT) 95% CI: [{ci_low:.3f}, {ci_high:.3f}], p={p_val:.4f}")


#percentage of senescenct cells according to hUSI and mayo
import pandas as pd
subset = adata.obs[adata.obs["comparison"].isin(["promet_AC", "nonpro_AC","promet_NT","nonpro_NT"])].copy()
subset["senescent"] = (subset["hUSI"] >= 0.6).astype(int)  # 1 = senescent, 0 = not
percentages = (
    subset.groupby("comparison")["senescent"]
    .value_counts(normalize=True)
    .rename("fraction")
    .reset_index()
)

percentages["percentage"] = percentages["fraction"] * 100
print(percentages)



#DEGS mayo 
with open(os.path.join(path_data,'AC_NT.pickle'), 'rb') as p:
    ac_nt = pickle.load(p)

with open(os.path.join(path_data,'pro_AC_vs_NT.pickle'), 'rb') as p:
    proac_nt = pickle.load(p)


dfs = []
categories= ['promet_AC_vs_NT']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = proac_nt.results[key]['df'].copy()
    df['label'] = cat
    dfs.append(df)

all_degs= pd.concat(dfs)
genes_of_interest = set(mayo_present) | {"PAEP"}   # union: all mayo genes + PAEP
df_mayo = all_degs.loc[all_degs.index.isin(genes_of_interest)].copy()


#volcano_plot
df_fixed= df_mayo.copy()
df_fixed['evidence'] = df_mayo['evidence'].replace(0, 1e-50)
df_fixed['evidence'] = df_fixed['evidence'].clip(lower=1e-50)
fig = volcano_plot_plot(
    df_fixed.query('comparison=="g0_vs_g1"'), effect_size='effect_size', evidence='evidence',
    n=30, annotate=True, xlim=(-2.5, 2.5),pseudocount=1e-50, title = "Promet AC vs Promet NT (Mayo)",
)
fig.tight_layout()
fig.savefig(os.path.join(path_results,"volcano_all_prometAC_vs_prometNT_mayo.png"),dpi=300)
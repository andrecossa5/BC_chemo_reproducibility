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
path_adata = os.path.join("/Users/ieo7295/Desktop/tests/cell/senescence_paep/data/SCT")
path_results= os.path.join(path_main, "results", "MDA")

#Data
adata= sc.read_h5ad(os.path.join(path_adata, "clustered.h5ad"))
adata_sc=sc.read(os.path.join(path_data, "clustered_norm.h5ad"))
sign_1 = pd.read_csv(os.path.join(path_data, "REACTOME_CELLULAR_SENESCENCE.v2024.1.Hs.tsv"), sep='\t')
sign_2 = pd.read_csv(os.path.join(path_data, "SAUL_SEN_MAYO.v2024.1.Hs.tsv"), sep='\t')
sign_3 = pd.read_csv(os.path.join("/Users/ieo7295/Downloads/FRIDMAN_SENESCENCE_UP.v2024.1.Hs.tsv"),sep='\t')

with open(os.path.join(path_data,'signatures.pickle'), 'rb') as f:
    data = pickle.load(f)


reactome = sign_1['REACTOME_CELLULAR_SENESCENCE'][16].split(',')
mayo = sign_2['SAUL_SEN_MAYO'][16].split(',')
friedman = sign_3['FRIDMAN_SENESCENCE_UP'][16].split(',')

reactome_present = []
for gene in reactome:
    if gene in adata.var_names:
        reactome_present.append(gene)

mayo_present = []
for gene in mayo:
    if gene in adata.var_names:
        mayo_present.append(gene)

friedman_present = []
for gene in friedman:
    if gene in adata.var_names:
        friedman_present.append(gene)

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


#add_column comparison
adata.obs['comparison'] = adata_sc.obs['comparison']
adata.uns['lognorm']
adata.X.max()

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
#     'fridman' : friedman_present,
#     'reactome' : reactome_present,
#     'sasp_sign' : sasp_sign,
#     'sasp_unique': sasp_unique

# }

# for key in d_sig:
#     gene_set = d_sig[key]
#     score_name = f"{key}_score"
#     gene_check = [g for g in gene_set if g in adata.var_names]
#     scores = scanpy_score(adata,gene_check)
#     adata.obs[f'{key}_score'] = scores.values

#cellula_scored signature.py
signatures = data['scores']
adata.obs['fridman_present'] = signatures['friedman_present']
adata.obs['mayo_present'] = signatures['mayo_present']
adata.obs['reactome_present'] = signatures['reactome_present']
adata.obs['PAEP'] = adata[:,'PAEP'].X.toarray().flatten()
adata.write_h5ad(os.path.join(path_adata, 'clustered_new.h5ad'))
#paep-senescence correlation
paep_expr = adata[:,'PAEP'].X.toarray().flatten()

df = pd.DataFrame({
    'PAEP': paep_expr,
    'mayo_present' : adata.obs['mayo_present'],
    'fridman_present' : adata.obs['fridman_present'],
    'reactome_present': adata.obs['reactome_present'], #    'sasp_sign_score': adata.obs['sasp_sign_score'],'sasp:unique_score': adata.obs['sasp_unique_score']
    'comparison': adata.obs['comparison']
})
colors = {
    'promet_AC': 'blue',
    'promet_NT': 'orange',
}

signatures = ['mayo_present','fridman_present','reactome_present','sasp:unique_score'] #'sasp_sign_score'
df = df[df['comparison'].isin(['promet_AC','promet_NT'])]
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
df['mayo_present'] = adata.obs['mayo_present']
df['UMAP1'] = adata.obsm['X_umap'][:, 0]
df['UMAP2'] = adata.obsm['X_umap'][:, 1]

sc.pl.umap(
    adata,
    color=['mayo_present', 'PAEP'],
    cmap='viridis',  
    vmax=[0.20, 6],
    title='Mayo Score and PAEP Co-localization'
)



import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 2, figsize=(10, 5))

# Mayo score with vmax
sc.pl.umap(
    adata,
    color='fridman_present',
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
plt.savefig(os.path.join(path_results, 'umap_mayo_paep'), dpi=300)



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

    if len(mayo) >= 30:  # avoid meaningless correlations
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

import scanpy as sc
import matplotlib.colors as clr
hUSI = cal_hUSI(adata)
adata.obs['hUSI'] = hUSI
sc.pp.highly_variable_genes(adata,n_top_genes=2000)
adata.raw = adata.copy()
adata = adata[:,adata.var.highly_variable]
sc.pp.pca(adata,n_comps=15)
sc.pp.neighbors(adata)
sc.tl.tsne(adata)
color_self = clr.LinearSegmentedColormap.from_list('pink_grey', ['#3AB370',"#EAE7CC","#FD1593"], N=256)
sc.pl.tsne(adata,color='hUSI',save='Python_demo_hUSI.png',size=30,cmap = color_self)


SenClass = SSE_hUSI(hUSI)

SenClass = GMM_hUSI(hUSI)

"""
Boxplots of DEGs from GSEA analysis.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib as plt
from typing import Dict, Iterable, Any, Tuple
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import plotting_utils as plu
from Cellula.plotting._plotting import *
from Cellula._utils import sanitize_neighbors

##


# Paths
path_main = '/Users/ieo7295/Desktop/BC_chemo_reproducibility/'
path_adata = os.path.join(path_main, 'data', 'MDA')
path_data = os.path.join(path_main, 'results', 'MDA','boxplot')
path_results = os.path.join(path_main, 'results', 'MDA', 'boxplot','plots')

#Data 
gl = pd.read_csv(os.path.join(path_data, "glycolisis_genes.txt"))
hp = pd.read_csv(os.path.join(path_data,"hypo_135.txt"))
mt = pd.read_csv(os.path.join(path_data,"mito_upregulated.txt"),sep='\t')
inf= pd.read_csv(os.path.join(path_data,"inflamm_29.txt"))
inf_ds= pd.read_csv(os.path.join(path_data,"inf_115.txt"))
cc = pd.read_csv(os.path.join(path_data,"cell_cycle_289.txt"))
dna_d = pd.read_csv(os.path.join(path_data,"DNA_damage.txt"))
inf_ss = pd.read_csv(os.path.join(path_data,"inflammatory_stress.txt"))
met_ss = pd.read_csv(os.path.join(path_data,"Metabolic_stress.txt"))
mit_ss = pd.read_csv(os.path.join(path_data,"Mitotic_stress.txt"))
promet_vs_nonpromet_NT = pd.read_csv(os.path.join(path_data,"promet_vs_nonpromet_NT_formatted.csv"))
promet_vs_nonpromet_AC = pd.read_csv(os.path.join(path_data,"promet_vs_nonpromet_AC_formatted.csv"))
adata= sc.read_h5ad(os.path.join(path_adata, "clustered.h5ad"))
adata_sc=sc.read(os.path.join(path_adata, "clustered_norm.h5ad"))
adata_ctc= sc.read(os.path.join("/Users/ieo7295/Desktop/tests/cell/data/default/clustered.h5ad"))

#add condition column to adata
adata.obs['comparison'] = adata_sc.obs['comparison']

#Viz

#1 boxplot of 56 glycolisis genes in prometastatic vs non-prometastatic clones (NT vs AC) 
glyco_genes= gl['Gene'].tolist()
gl_gene= [x for x in glyco_genes if x in adata.var_names]
gl_gene = list(set(gl_gene))
adata_gl = adata[:,gl_gene].copy()

adata_gl.obs['mean_expr'] = adata_gl.X.mean(axis=1).A1 

df = adata_gl.obs[['mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample', 'comparison'])['mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(5.3,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','promet_AC'],
           ['nonpro_NT','nonpro_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 52 glycolysis genes', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_52_glycolisis_mean_expr.png'), dpi=300)


#2 boxplot : 13 hypoxia common and 134 non-redundant in three condition (promet_NT, nonpro_AC, promet_AC)
hypo_13_genes = ["ANGPTL4","BTG1","DUSP1","EFNA1","ERRFI1","ETS1","IL6","KDM3A"
                 ,"MT2A","P4HA1","PLIN2","TNFAIP3","TPBG"]

hypo_gene_13= [x for x in hypo_13_genes if x in adata.var_names]
hypo_gene_13 = list(set(hypo_gene_13))
adata_hypo_13 = adata[:,hypo_gene_13].copy()

adata_hypo_13.obs['hypo_mean_expr'] = adata_hypo_13.X.mean(axis=1).A1 

df = adata_hypo_13.obs[['hypo_mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample', 'comparison'])['hypo_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(5.3,5))
order = ['promet_NT','nonpro_AC','promet_AC']

plu.box(df_bulk,
    x='comparison', y='hypo_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='hypo_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 13 hypoxia genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_13_hypoxia_sample.png'), dpi=300)

#134
hp_genes = hp['Genes'].tolist()

hp_gene_134= [x for x in hp_genes if x in adata.var_names]
hp_gene_134 = list(set(hp_gene_134))
adata_hp_134 = adata[:, hp_gene_134].copy()

adata_hp_134.obs['hypo_mean_expr'] = adata_hp_134.X.mean(axis=1).A1 

df = adata_hp_134.obs[['hypo_mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample', 'comparison'])['hypo_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(5.3,5))
order = ['promet_NT','nonpro_AC','promet_AC']

plu.box(df_bulk,
    x='comparison', y='hypo_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='hypo_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 134 hypoxia genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_134_hypoxia_sample.png'), dpi=300)


#mito upregulated (nonpro_NT,nonpro_AC,promet_NT,promet_AC) , 124 common
genes_124 = mt[
    (mt['Treated Prometastatic'] == 'Yes') & 
    (mt['Treated NonPrometastastic'] == 'Yes')
]['Gene'].tolist()

mt_gene_124= [x for x in genes_124 if x in adata.var_names]
mt_gene_124 = list(set(mt_gene_124))
adata_mt_124 = adata[:, mt_gene_124]

adata_mt_124.obs['mito_mean_expr'] = adata_mt_124.X.mean(axis=1).A1 

df = adata_mt_124.obs[['mito_mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample', 'comparison'])['mito_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(6.9,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mito_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mito_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expr of 124 common (promet_AC & nonpro_AC) mitochondrial genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_124_mitochondrial_sample.png'), dpi=300)

#241 promet_AC
genes_220 = mt[
    (mt['Treated Prometastatic'] == 'Yes')
]['Gene'].tolist()

mt_gene_220= [x for x in genes_220 if x in adata.var_names]
mt_gene_220 = list(set(mt_gene_220))
mt_gene_220 = adata[:, mt_gene_220].copy()

mt_gene_220.obs['mito_mean_expr'] = mt_gene_220.X.mean(axis=1).A1 

df = mt_gene_220.obs[['mito_mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample', 'comparison'])['mito_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(6.9,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mito_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mito_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expr of 220 promet_AC mitochondrial genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_220_mitochondrial_sample.png'), dpi=300)

#Inflammation 29 promet_NT 
genes_29_inf = inf['Genes'].tolist()

inf_gene_29= [x for x in genes_29_inf if x in adata.var_names]
inf_gene_29 = list(set(inf_gene_29))
adata_inf_29 = adata[:, inf_gene_29]

adata_inf_29.obs['inf_mean_expr'] = adata_inf_29.X.mean(axis=1).A1 

df = adata_inf_29.obs[['inf_mean_expr', 'sample', 'comparison','condition']].copy()
df_bulk = df.groupby(['sample','comparison'])['inf_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(6.9,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='inf_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='inf_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 29 inflammation genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_29_inf_sample.png'), dpi=300)

#Inflammation 115 common (promet_AC & nonpro_AC)
genes_115_inf = inf_ds['Genes'].tolist()

inf_gene_114= [x for x in genes_115_inf if x in adata.var_names]
inf_gene_114= list(set(inf_gene_114))
adata_inf_114 = adata[:, inf_gene_114]

adata_inf_114.obs['inf_mean_expr'] = adata_inf_114.X.mean(axis=1).A1 

df = adata_inf_114.obs[['inf_mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['inf_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='inf_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='inf_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 114 common inflammation genes (promet_AC & nonpro_AC)', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_114_inf_sample.png'), dpi=300)




#stress_response: DNA damage
dna_d_genes = dna_d['Gene'].tolist()

dna_d_50= [x for x in dna_d_genes if x in adata.var_names]
dna_d_50= list(set(dna_d_50))
adata_dna_d_50 = adata[:, dna_d_50]

adata_dna_d_50.obs['dna_mean_expr'] = adata_dna_d_50.X.mean(axis=1).A1 

df = adata_dna_d_50.obs[['dna_mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['dna_mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='dna_mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='dna_mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 50 DNA damage genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_50_dnad_sample.png'), dpi=300)

#stress response: ER manca ALR
er_genes= ["ERO1A","OSTC","ERP29","QSOX1","RPN1","SDF2L1","SEC11A",
            "SERP1","ALR","EDEM1","HSP90B1","SVIP","UFM1"]
len(er_genes)
er_514= [x for x in er_genes if x in adata.var_names]
er_514= list(set(er_514))
adata_er_514 = adata[:, er_514]

adata_er_514.obs['mean_expr'] = adata_er_514.X.mean(axis=1).A1 

df = adata_dna_d_50.obs[['mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 50 DNA damage genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_50_dnad_sample.png'), dpi=300)

#hypoxia stress
hy_gene= ["ANGPTL4","EGLN2","FAM162A","NRP1","P4HA1","PGF",
            "PLOD2","VEGFA","CCN2","ELK3","ENG","MALAT1",
            "LDHA","PGK1","BHLHE40","DDIT4"]

hy_16= [x for x in hy_gene if x in adata.var_names]
hy_16= list(set(hy_16))
adata_hy_16 = adata[:, hy_16].copy()

adata_hy_16.obs['mean_expr'] = adata_hy_16.X.mean(axis=1).A1 

df = adata_hy_16.obs[['mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 16 Hypoxia stress genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_16_hypoxia_stress_sample.png'), dpi=300)

#stress : Inflammatory 
inf_ss_genes = inf_ss['Gene'].tolist()

inf_ss_32= [x for x in inf_ss_genes if x in adata.var_names]
inf_ss_32= list(set(inf_ss_32))
adata_inf_ss_32 = adata[:, inf_ss_32].copy()

adata_inf_ss_32.obs['mean_expr'] = adata_inf_ss_32.X.mean(axis=1).A1 

df = adata_inf_ss_32.obs[['mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 32 Inflammation stress genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_32_inflammation_stress_sample.png'), dpi=300)

#Stress: metabolic 37 
met_ss_genes = met_ss['Gene'].tolist()

met_ss_38= [x for x in met_ss_genes if x in adata.var_names]
met_ss_38= list(set(met_ss_38))
adata_met_ss_38 = adata[:, met_ss_38].copy()

adata_met_ss_38.obs['mean_expr'] = adata_met_ss_38.X.mean(axis=1).A1 

df = adata_met_ss_38.obs[['mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','nonpro_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 32 Inflammation stress genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_32_inflammation_stress_sample.png'), dpi=300)

#stress: mitotic
mit_ss_genes = mit_ss['Gene'].tolist()

mit_ss_24= [x for x in mit_ss_genes if x in adata.var_names]
mit_ss_24= list(set(mit_ss_24))
adata_mit_ss_24 = adata[:, mit_ss_24].copy()

adata_mit_ss_24.obs['mean_expr'] = adata_mit_ss_24.X.mean(axis=1).A1 

df = adata_mit_ss_24.obs[['mean_expr', 'sample', 'comparison']].copy()
df_bulk = df.groupby(['sample','comparison'])['mean_expr'].mean().reset_index()

fig, ax = plt.subplots(figsize=(7.2,5))
order = ['nonpro_NT','nonpro_AC','promet_NT','promet_AC']

plu.box(df_bulk,
    x='comparison', y='mean_expr', ax=ax, color='grey', add_stats=True, 
    pairs=[['promet_NT','nonpro_AC'],['promet_NT','promet_AC'],['nonpro_AC','promet_AC'],['nonpro_NT','promet_NT'],['nonpro_NT','promet_AC'],
           ['nonpro_NT','nonpro_AC']], 
    x_order=order
)
plu.strip(df_bulk, x='comparison', y='mean_expr', ax=ax, color='k', x_order=order)
plu.format_ax(ax=ax, title='Mean expression of 24 Mitotic stress genes', ylabel='mean', rotx=90, reduced_spines=True)
fig.tight_layout()
fig.savefig(os.path.join(path_results, f'Boxplot_24_mitotic_stress_sample.png'), dpi=300)








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

#cell_cycle 
cc_genes = cc['Genes'].tolist()
adata_cc = adata[adata.obs['comparison'].isin(['promet_NT','promet_AC'])].copy()
gene_list = [x for x in cc_genes if x in adata_cc.var_names]


sc.tl.score_genes(adata_cc, gene_list=gene_list, score_name='prolif_score')

#prolif_score
fig, ax = plt.subplots(figsize=(7.2,5))
order= ['promet_NT','promet_AC']
pairs= [('promet_NT','promet_AC')]
violin(adata_cc.obs, 'comparison', 'prolif_score', ax=ax, color='darkgrey', add_stats=True, x_order=order,pairs=pairs)
plu.format_ax(ax, title='Proliferation Score in prometastatic cells', ylabel='Score')
plt.tight_layout()
fig.savefig(os.path.join(path_results, 'cc_prolif_promet.png'), dpi=400)

#cycling
fig = plt.figure(figsize=(14,7))
gs = GridSpec(2, 6, figure=fig, height_ratios=[1,1.5])
ax = fig.add_subplot(gs[0,1:-1])
pairs = [('promet_NT','promet_AC')]

embs = (
    adata_cc.obs
    .join(
    pd.DataFrame(
        adata_cc.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata_cc.obs_names
    ))
)

order= ['promet_NT','promet_AC']
violin(adata_cc.obs, 'comparison', 'cycling', ax=ax, color='darkgrey', add_stats=True, x_order=order,pairs=pairs)
plu.format_ax(ax, title='Cell cycle signatures scores', ylabel='Score')
ax.spines[['left', 'right', 'top']].set_visible(False)

ax = fig.add_subplot(gs[1,:2])
draw_embeddings(embs, cont='s_seurat', ax=ax, title='cycling', cbar_kwargs={'palette':'mako'}, s=1)
ax.axis('off')
ax = fig.add_subplot(gs[1,2:4])
draw_embeddings(embs, cont='G1/S', ax=ax, title='G1/S', cbar_kwargs={'palette':'mako'}, s=1)
ax.axis('off')
ax = fig.add_subplot(gs[1,4:])
draw_embeddings(embs, cont='G2/M', ax=ax, title='G2/M', cbar_kwargs={'palette':'mako'}, s=1)
ax.axis('off')

fig.tight_layout()
fig.savefig(os.path.join(path_results, 'cc_promet.png'), dpi=400)



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
from scipy.sparse import csr_matrix
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *
from Cellula.plotting._plotting import *
from matplotlib.gridspec import GridSpec
from plotting_utils._plotting_base import *
from Cellula.plotting._plotting import *
from Cellula._utils import sanitize_neighbors
from BC_chemo_utils.utils import prep_df_for_dotplot
from BC_chemo_utils.tests import *
from BC_chemo_utils.plotting import *
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


## DE ##

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata_auc, path_data, contrasts_name='paep_contrasts')

# Here we go
D = Dist_features(adata_auc, contrasts, jobs=jobs, app=False)
D.run_all_jobs()

dfs = []
categories= ['nonpro_AC_vs_NT', 'promet_AC_vs_NT', 'promet_AC', 'pro_nonpro_AC', 'promet_NT']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, "Degs_regulon.csv"))


#Embeddings regulons

# UMAP #15 0,4
runUmap = umap.UMAP(n_neighbors=30, min_dist=0.4, metric='euclidean').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv(os.path.join(path_results, "scenic_umap_30.txt"), sep='\t')

dr_umap = pd.read_csv(os.path.join(path_results,'scenic_umap_30.txt'), sep='\t', header=0, index_col=0 )

#Viz 
umap_array = dr_umap[['X', 'Y']].values
adata_auc.obsm['X_umap'] = umap_array

embs = (
    adata_auc.obs
    .join(
    pd.DataFrame(
        adata_auc.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata_auc.obs_names
    ))
)
  
df_markers = pd.read_csv(os.path.join(path_data, 'Degs_regulon.csv'), index_col=0)



# Create colors
cell_state_colors = create_palette(embs, 'condition', 'tab20')
cell_state_colors['Undefined'] = '#E8E7E7'

## UMAP cell_state + dotplot


#def prep_df_for_dotplot(df_markers, n=3):
"""
Prep df for dotplot.
"""


genes = {}
comparisons = df_markers['comparison'].unique()
n=3
for comparison in comparisons:
    group = '_vs_'.join(comparison.split('_vs_')[:2])                         
    print(group)
    genes[group] = (
        df_markers.query('comparison == @comparison and perc_FC>1 and AUROC>0.8')
        .index[:n]
        .to_list()
    )
    print(genes)
from itertools import chain
genes_ = list(chain.from_iterable([ genes[x] for x in genes ]))

# Get percs and log2FC
df = (
    df_markers
    .loc[genes_, ['effect_size', 'group_perc', 'comparison']]
    .reset_index()
    .rename(columns={'index':'regulon', 'effect_size':'log2FC'})
)

#Clip 
df['log2FC'][df['log2FC'] <= np.percentile(df['log2FC'], 5)] = np.percentile(df['log2FC'], 5)
df['log2FC'][df['log2FC'] >= np.percentile(df['log2FC'], 95)] = np.percentile(df['log2FC'], 95)

    #return df






fig = plt.figure(figsize=(20,10))
gs = GridSpec(1, 2, figure=fig, width_ratios=[2,2.5])

ax = fig.add_subplot(gs[0])
draw_embeddings(
    embs, cat='condition', ax=ax, title='Condition', 
    legend_kwargs={
        'ncols':1, 'colors':cell_state_colors, 
        'bbox_to_anchor':(1,1), 'loc':'upper left'
    },
)

ax = fig.add_subplot(gs[1])
df_ = prep_df_for_dotplot(df_markers)
df['comparison'] = df['comparison'].map(lambda x: x.join('_')[0][:-1]) 

sns.scatterplot(
    data=df, y='comparison', x='regulon', size='group_perc', hue='log2FC', 
    palette='mako', ax=ax, sizes=(1, 100)
)
format_ax(ax, title='Markers', xlabel='Top 3 marker regulons', ylabel='condition', rotx=90)
ax.legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

fig.subplots_adjust(left=.05, right=.9, top=.9, bottom=.2, wspace=1.3)
# Save
fig.savefig(os.path.join(path_results, 'regulon_condition_30.png'), dpi=500)










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
        if tf is None:
            print("Missing 'regulon' key in item:", item)
            continue

        if tf not in alt_mtx.columns:
            print(f"Warning: {tf} column is not found in alt_mtx.")
            continue

        motif_data = item.get('motifData', {})
        print(f"Motif Data: {motif_data}")

        genes_for_set = alt_mtx.index[alt_mtx[tf] == 1].tolist()
        print(f"Genes for {tf}: {genes_for_set}")

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
reg_paep = {}
for tf, value in d_reg.items():
    gene_set = value.get('gene_set', [])
    if 'PAEP' in gene_set:
        reg_paep[tf] = {
            'gene_set': gene_set,
            'motif_data': value.get('motif_data', {})
        }


#Viz 
jobs,contrasts= prep_jobs_contrasts(adata_auc, path_data, contrasts_name='paep_contrasts')

dist = Dist_features(adata_auc, contrasts=contrasts, jobs=jobs)
dist.select_genes()  

#create column comparison
contrasts
contrast= dist.contrasts['pro_nonpro_AC']

adata_auc.obs['pro_nonpro_AC'] = contrast.category

cols_to_merge = ['nonpro_AC_vs_NT', 'promet_AC_vs_NT', 'promet_AC', 'promet_NT']

# Apply row-wise logic to keep the non-'to_exclude' value
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

adata_auc_thrb_clean = adata_auc_thrb[adata_auc_thrb.obs['comparison'].notna()].copy()


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
#ax.spines[['left', 'right', 'top']].set_visible(False)
fig.tight_layout()

fig.savefig(os.path.join(path_results, 'THRB_violin.png'), dpi=400)




sns.boxplot(x='comparison', y='THRB', data=adata_auc_thrb_clean.obs)
plt.show()
sns.violinplot(x='comparison', y='THRB', data=adata_auc_thrb_clean.obs)
plt.show()

adata_auc_thrb.var.drop('THRB', inplace=True)

#umap
adata_tmp=adata.copy()
adata_tmp.obs['TF_THRB']=adata_auc_thrb.obs['THRB']
adata_tmp.obs['comparison']= adata_auc_thrb.obs['comparison']
regulon_to_plot=['TF_THRB']

#save
plt.figure(figsize=(8, 6))  
sc.pl.umap(adata_tmp, color='comparison', cmap='viridis', vmin=0, show=False)
plt.savefig("umap_comparison.png", dpi=300, bbox_inches='tight')  
plt.close()

plt.figure(figsize=(8, 6))  
sc.pl.umap(adata_tmp, color=regulon_to_plot, cmap='viridis', vmin=0, show=False)
plt.savefig("umap_THRB.png", dpi=300, bbox_inches='tight')  
plt.close()





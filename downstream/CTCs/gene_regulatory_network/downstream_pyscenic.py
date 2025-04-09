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
adata_auc.var['regulon'] = auc_mtx.columns
adata_auc.var['regulon'] = adata_auc.var['regulon'].str.replace('(+)', '', regex=False)
adata_auc.var_names = adata_auc.var_names.str.replace('(+)', '',regex=False)

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
categories= ['nonpro_AC_vs_NT', 'promet_AC_vs_NT', 'promet_AC', 'dom_promet_vs_nonmet', 'promet_NT']

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, "Degs_regulon.csv"))



#Embeddings regulons

# UMAP
# runUmap = umap.UMAP(n_neighbors=15, min_dist=0.4, metric='correlation').fit_transform
# dr_umap = runUmap( auc_mtx )
# pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv(os.path.join(path_results, "scenic_umap.txt"), sep='\t')

dr_umap = pd.read_csv(os.path.join(path_results,'scenic_umap.txt'), sep='\t', header=0, index_col=0 )

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
    group = '_'.join(comparison.split('_')[:4])
    print(group)
    genes[group] = (
        df_markers.query('comparison == @comparison and perc_FC>1 and AUROC>0.8')
        .index[:n]
        .to_list()
    )

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










fig = plt.figure(figsize=(20,8))
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
df['comparison'] = df['comparison'].map(lambda x: x.join('_')[3][:-1]) 

sns.scatterplot(
    data=df_, y='comparison', x='gene', size='group_perc', hue='log2FC', 
    palette='mako', ax=ax, sizes=(1, 100)
)
format_ax(ax, title='Markers', xlabel='Top 3 marker regulons', ylabel='condition', rotx=90)
ax.legend(loc='center left', bbox_to_anchor=(1,.5), frameon=False)

fig.subplots_adjust(left=.05, right=.9, top=.9, bottom=.2, wspace=1.3)
plt.show()
# Save
fig.savefig(os.path.join(path_results, 'annotated_cell_states.png'), dpi=500)










#Reformat regulons 
regulons.columns = regulons.iloc[0]         
regulons = regulons.drop(index=[0])     
regulons.columns.values[0] = regulons.iloc[0, 0]  
regulons.columns.values[1] = regulons.iloc[0, 1]
regulons= regulons.drop(index=[1]).reset_index(drop=True)

#Regulon : tf + gene_set
alt_mtx = pd.DataFrame(lf.ra.Regulons, index=lf.ra.Gene)


#Create dict
d_reg={}
for tf in alt_mtx.columns:
    genes_for_set = alt_mtx.index[alt_mtx[tf] == 1].tolist()
    d_reg[tf] = genes_for_set


reg_paep={}
#N regulons with PAEP
for key, value in d_reg.items():
    for val in value: 
        if val == 'PAEP':
            reg_paep[key]= value


#
df_motif
rank_df_10
rank_df_500
regulons





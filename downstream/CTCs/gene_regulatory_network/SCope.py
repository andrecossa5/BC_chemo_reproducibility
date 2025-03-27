"""
Create loom file for SCope visualization
"""

import os
import json
import zlib
import base64
import umap
import pandas as pd
import loompy 
from openTSNE import TSNE




#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility" 
path_data = os.path.join(path_main,"data", "CTCs", "grn", "pyscenic_res") 
path_results = os.path.join(path_main,"data", "CTCs", "grn", "results")

#Data
auc_mtx = pd.read_csv(os.path.join(path_data, "ctc_aucell.csv"),index_col=0)
X = auc_mtx.to_numpy() 


# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv(os.path.join(path_results, "scenic_umap.txt"), sep='\t')
# tSNE
tsne = TSNE(perplexity=30, n_jobs=-1)
X_tsne = tsne.fit(X)
pd.DataFrame(X_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv(os.path.join(path_results, "scenic_tsne.txt"), sep='\t')

#Here we go

row_attrs = {"Regulons": auc_mtx.index.values}  
col_attrs = {"CellID": auc_mtx.columns.values} 

loom_file = os.path.join(path_results,"auc_mtx.loom")

loompy.create(loom_file, auc_mtx.values, row_attrs, col_attrs)

loom_file = os.path.join(path_results,"auc_mtx.loom")
ds = loompy.connect(loom_file, mode="r") 
dr_umap = pd.read_csv(os.path.join(path_results, 'scenic_umap.txt'))
dr_tsne = pd.read_csv(os.path.join(path_results, 'scenic_tsne.txt'))


scope_loom = export_to_loom.SCopeLoom.read_loom(filename=args.loom_input)
scope_loom.add_embedding(embedding=dr_umap, embedding_name="SCENIC AUC UMAP", is_default=True)
scope_loom.add_embedding(embedding=dr_tsne, embedding_name="SCENIC AUC t-SNE", is_default=False)
#scope_loom.export(out_fname=)

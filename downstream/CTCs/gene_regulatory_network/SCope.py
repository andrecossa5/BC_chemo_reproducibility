"""
Create loom file for SCope visualization
"""

import os
import json
import zlib
import base64
import umap
import pandas as pd
#from sklearn.manifold import TSNE
from openTSNE import TSNE


#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility" 
path_data = os.path.join(path_main,"data", "CTCs", "grn", "pyscenic_res") 
path_results = os.path.join(path_main,"data", "CTCs", "grn", "results")

#Data
auc_mtx = pd.read_csv(os.path.join(path_data, "ctc_aucell.csv"),index_col=0)
regulons = pd.read_csv(os.path.join(path_data, 'ctc_regulons.csv'))

#UMAP 
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv(os.path.join(path_results, "scenic_umap.txt"),sep='\t')
# tSNE
X = auc_mtx.to_numpy()
tsne = TSNE(
    perplexity=30, 
    n_jobs=-1  # Use all CPU cores
)
X_tsne = tsne.fit(X)
pd.DataFrame(X_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv(os.path.join(path_results, "scenic_tsne.txt"),sep='\t')


dr_umap = pd.read_csv(os.path.join(path_results, 'scenic_umap.txt'), sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv(os.path.join(path_results, 'scenic_tsne.txt'), sep='\t', header=0, index_col=0 )

#Fix regulon objects to display properly in SCope
auc_mtx.columns = auc_mtx.columns.str.replace("\(",'_(')
regulons.columns = [x.replace("(","_(") for x in regulons.columns]


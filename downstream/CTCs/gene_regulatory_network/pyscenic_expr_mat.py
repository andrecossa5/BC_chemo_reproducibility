"""
Create loom file for pyscenic (raw)
"""

import os 
import pandas as pd
import numpy as np
import loompy as lp
import scanpy as sc

#Files required by pyscenic 

# "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
# "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
# "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
# "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather"
# "allTFs_hg38.txt"
# "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"


#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility"
path_data= os.path.join(path_main,"data", "CTCs")
path_tfs=os.path.join(path_main,"data", "CTCs", "grn", "resources", "allTFs_hg38.txt")


#output path
f_loom_path_scenic= os.path.join(path_main,"data", "CTCs", "grn", "resources", "ctc_scenic.loom")

#Data
adata= sc.read(os.path.join(path_data, "clustered.h5ad" ))
tfs= [tf.strip() for tf in open(path_tfs)]


# as a general QC. We inspect that our object has transcription factors listed in our main annotations.
print(
    f"%{np.sum(adata.var.index.isin(tfs))} out of {len(tfs)} TFs are found in the object"
)

#Select raw counts
adata.X=adata.layers['raw']

# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)

loom_file= os.path.join(path_data, 'grn', 'resources', 'ctc_scenic_raw.loom')
ds= lp.connect(loom_file,mode='r')




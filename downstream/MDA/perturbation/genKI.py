#!/usr/bin/python

"""
Perturbation analysis.
"""

import os
import argparse
import scanpy as sc
import GenKI as gk
from scipy.sparse import csr_matrix
from GenKI.dataLoader import DataLoader
from GenKI.train import VGAE_trainer
from GenKI import utils
from plotting_utils._plotting_base import *
from gseapy import enrichr
matplotlib.use('macOSX')


##


my_parser = argparse.ArgumentParser(
    prog='genKI',
    description=
    """
    Gene KI script.
    """
)

# Input
my_parser.add_argument(
    '--path_main', 
    type=str,
    default=None,
    help='Path_main. Default: None.'
)

# path_GRN
my_parser.add_argument(
    '--path_GRN', 
    type=str,
    default=None,
    help='Path_main. Default: None.'
)

# Output
my_parser.add_argument(
    '--gene', 
    type=str,
    default=None,
    help='Gene to KO. Default: None.'
)

# treshold
my_parser.add_argument(
    '--n_ORA',
    type=int,
    default=50,
    help='n gene for ORA. Default: None.'
)

# treshold
my_parser.add_argument(
    '--gene_set',
    type=str,
    default='MSigDB_Hallmark_2020',
    help='Gene set for ORA. Default: None.'
)

# Markers
my_parser.add_argument( 
    '--rebuild_GRN', 
    action='store_true',
    help='Rebuild GRN. Default: True.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_main = args.path_main
path_GRN = args.path_GRN
gene = args.gene
n_ORA = args.n_ORA
gene_set = args.gene_set


##


# Paths
# path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'perturbation')

# # Define params
# gene = 'PAEP'
# rebuild_GRN = False
# gene_set = 'MSigDB_Hallmark_2020'
# n = 50


##


def main():

    # Read clustered adata
    adata = sc.read(os.path.join(path_data, 'subsample.h5ad'))

    # Select HVGs, and scale 
    genes = adata.var['highly_variable_features']
    adata = adata[:,genes].copy()
    adata.layers['norm'] = csr_matrix(adata.X)
    sc.pp.scale(adata)
    adata.X = csr_matrix(adata.X)

    # GRN construction
    data_wrapper = DataLoader(adata, target_gene=[gene], rebuild_GRN=args.rebuild_GRN, 
                            GRN_file_dir=path_GRN, pcNet_name="pcNet", n_cpus=8)
    data_wt = data_wrapper.load_data()
    data_ko = data_wrapper.load_kodata()

    # Train
    hyperparams = {"epochs": 100, "lr": 7e-4,  "beta": 1e-4, "seed": 1234}
    sensei = VGAE_trainer(
        data_wt, epochs=hyperparams["epochs"], lr=hyperparams["lr"], 
        log_dir=None, beta=hyperparams["beta"], seed=hyperparams["seed"]
    )
    sensei.train()

    # Rank KO-affected genes
    z_mu_wt, z_std_wt = sensei.get_latent_vars(data_wt)
    z_mu_ko, z_std_ko = sensei.get_latent_vars(data_ko)
    dis = gk.utils.get_distance(z_mu_ko, z_std_ko, z_mu_wt, z_std_wt, by="KL")
    res_raw = utils.get_generank(data_wt, dis, rank=True)

    # Get significantly affected by permutation test
    null = sensei.pmt(data_ko, n=100, by="KL")
    gene_list_df = utils.get_generank(data_wt, dis, null)

    # ORA
    ORA_df = enrichr(
        gene_list=gene_list_df.index.to_list()[:n_ORA],
        gene_sets=[gene_set],
        organism='human', 
        outdir=None, 
    ).results

    # Save
    os.mkdir(os.path.join(path_results, gene))
    gene_list_df.to_csv(os.path.join(path_results, gene, 'KO_genes.csv'))
    ORA_df.set_index('Gene_set').to_csv(os.path.join(path_results, 'KO_genes_ORA.csv'), index_col=0)


##


# GRN loading
# GRN = np.load('GRNs/pcNet.npz')
# X = csr_matrix((GRN['data'], GRN['indices'], GRN['indptr']), shape=GRN['shape'])
# df = pd.DataFrame(X.A, index=genes, columns=genes)
# 
# fig, ax = plt.subplots(figsize=(6,6))
# f = lambda x: np.where(x>0.00001,1,0)
# L = linkage(f(df), method='weighted')
# order = leaves_list(L)
# ax.imshow(f(df)[np.ix_(order, order)])
# fig.tight_layout()
# plt.show()

# d = {}
# for n in range(50, 100, 10):
#     genes = df.median(axis=1).sort_values(ascending=False).index[:n]
#     d[n] = enrichr(
#         gene_list=genes.to_list(),
#         gene_sets=['GO_Biological_Process_2021'],
#         organism='human', 
#         outdir=None, 
#     ).results['Term'][:10]


# Run
if __name__ == "__main__":
    main()


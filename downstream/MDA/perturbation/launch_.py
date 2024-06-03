#!/usr/bin/python

import os
import pandas as pd

path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_GRN = os.path.join(path_main, 'results/perturbation/GRNs')
gene_list = pd.read_csv(os.path.join(path_main, 'data/MDA/GenKI_genes_200524.csv'), index_col=0).index

for gene in gene_list:
    print(f'Gene {gene}...')
    pycall = f'python genKI.py --path_main {path_main} --path_GRN {path_GRN} --gene {gene}'
    os.system(pycall)



"""
Script to perform GSEA and visualization of pseudobulk DE.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
from itertools import chain, product
from functools import reduce
from plotting_utils import *
from BC_chemo_utils.tests import fastGSEA
from BC_chemo_utils.utils import read_model_results, read_GSEA
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_results = os.path.join(path_main, 'results', 'MDA', 'pseudobulk')

# Read DE results
df = (
    read_model_results(path_results, 'fixed')
    [['logFC', 'FDR', 'contrast', 'regressed_cc']]
    .assign(model='fixed')
)


##


# GSEA on all DEGs lists
combos = product(
    df['contrast'].unique(), 
    df['model'].unique(), 
    df['regressed_cc'].unique()
)
for c, m, cc in combos:
    df_ = fastGSEA(df.query('contrast==@c and model==@m and regressed_cc==@cc')['logFC'])
    df_.to_csv(
        os.path.join(
            path_results, 'GSEA', 'GO_biological_process_2021', f'{c}_{m}_{cc}.csv'
        )
    )

# Read GSEA
df_gsea = read_GSEA(path_results)


##


# Make some sense out of it!

# Top represented genes across all contrasts: within top5, FDR and logFC
feat = 'FDR'
pd.Series(chain.from_iterable(
    df
    .groupby(['contrast', 'model', 'regressed_cc'])
    .apply(lambda x: x.sort_values(feat, ascending=True).index.to_list()[:5])
)).value_counts().to_frame('n').to_csv(
    os.path.join(path_results, f'top_represented_genes_{feat}.csv')
)

# Over-represented terms and genes in GSEA

# Terms
feat = 'p_adjusted'
pd.Series(chain.from_iterable(
    df_gsea
    .groupby(['contrast', 'model', 'cc'])
    .apply(lambda x: x.sort_values(feat, ascending=True).index.to_list()[:5])
)).value_counts().to_frame('n').to_csv(
    os.path.join(path_results, f'top_represented_terms_{feat}.csv')
)

# Leading genes
pd.Series(chain.from_iterable(
    df_gsea
    .groupby(['contrast', 'model', 'cc'])
    .apply(lambda x: x.head(5)['Lead_genes'].map(lambda x: x.split(';')))
)).value_counts().to_frame('n').to_csv(
    os.path.join(path_results, f'top_represented_leading_genes.csv')
)


##
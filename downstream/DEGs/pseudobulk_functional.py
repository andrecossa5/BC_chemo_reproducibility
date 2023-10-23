"""
Script to perform GSEA and visualization of pseudobulk DE.
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
from itertools import product
from plotting_utils import *
from BC_chemo_utils.tests import fastGSEA
matplotlib.use('macOSX')


#


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_results = os.path.join(path_main, 'results', 'MDA', 'pseudobulk')

# Read DE results
_fixed = []
for path_root, _, files in os.walk(os.path.join(path_results, 'fixed')):
    for x in files:
        _fixed.append(os.path.join(path_root, x))
_mixed = []
for path_root, _, files in os.walk(os.path.join(path_results, 'random')):
    for x in files:
        _mixed.append(os.path.join(path_root, x))

# Merge all outputs in a single-table
DE_fixed_l = []
for x in _fixed:
    DE_fixed_l.append(
        pd.read_csv(x, index_col=0)
        .assign(regressed_cc=bool(re.search('G2.M', x)))
    )
DE_fixed = pd.concat(DE_fixed_l)[['logFC', 'FDR', 'contrast', 'regressed_cc']].assign(model='fixed')

DE_mixed_l = []
for x in _fixed:
    DE_mixed_l.append(
        pd.read_csv(x, index_col=0)
        .assign(regressed_cc=bool(re.search('G2.M', x)))
    )
DE_mixed = (
    pd.concat(DE_mixed_l)
    .rename(columns={'adj.P.Val':'FDR'})
    [['logFC', 'FDR', 'contrast', 'regressed_cc']]
    .assign(model='mixed')
)

# Merge and save
df = pd.concat([ DE_fixed, DE_mixed ])
df.to_csv(os.path.join(path_results, 'results.csv'))


##


# GSEA on all DEGs lists
combos = list(product(df['contrast'].unique(), df['model'].unique(), df['regressed_cc'].unique()))
for c, m, cc in combos:
    df_ = fastGSEA(df.query('contrast==@c and model==@m and regressed_cc==@cc')['logFC'])
    df_.to_csv(
        os.path.join(path_results, 'GSEA', 'GO_biological_process_2021', '_'.join([c, m, str(cc), '.csv'])
    ))


##


# Make some sense out of it!

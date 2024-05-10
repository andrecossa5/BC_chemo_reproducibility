"""
Volcano and lollipop plots of pseudobulk clones.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import scanpy as sc
from matplotlib.gridspec import GridSpec
from plotting_utils._plotting_base import *
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'pseudobulk')


##


# Read DE and volcano plots
name_files = [
    'results_treated_vs_untreated_PTs_nUMIs_mito_perc_G2.M.csv',
    'results_double_vs_single_treated_lungs_nUMIs_mito_perc_G2.M.csv',
    'results_double_vs_untreated_lungs_nUMIs_mito_perc_G2.M.csv',
    'results_single_vs_untreated_lungs_nUMIs_mito_perc_G2.M.csv'
]
names = [ 
    'PT, treated vs untreated',
    'lung, double- vs single-treated',
    'lung, double- vs untreated',
    'lung, single- vs untreated'
]
d = {}
for name, name_file in zip(names, name_files):
    d[name] = pd.read_csv(os.path.join(path_results, 'fixed', name_file), index_col=0)

# Volcanoes
for k in d:
    fig = volcano(d[k], effect_size='logFC', evidence='FDR', n=5, annotate=True, xlim=(-5, 5))
    fig.tight_layout()
    fig.suptitle(k)
    fig.savefig(os.path.join(path_results, 'volcano', f'{k}.png'), dpi=300)


##


# Read GSEA and lollipops plots
path_ = os.path.join(path_results, 'GSEA', 'GO_biological_process_2021')
name_files = [
    'treated_vs_untreated_PTs_fixed_True.csv',
    'double_vs_single_treated_lungs_fixed_True.csv',
    'double_vs_untreated_lungs_fixed_True.csv',
    'single_vs_untreated_lungs_fixed_True.csv'
]
names = [ 
    'PT, treated vs untreated',
    'lung, double- vs single-treated',
    'lung, double- vs untreated',
    'lung, single- vs untreated'
]
d = { name : pd.read_csv(os.path.join(path_, filename)) for name,filename in zip(names, name_files) }

# Plot
for x in d:  
    fig, ax = plt.subplots(figsize=(6,5))
    stem_plot(d[x].head(10), 'NES', ax=ax)
    format_ax(ax, xlabel='NES', yticks=d[x].tail(10)['Term'].map(lambda x: x.split('(')[0] ))
    fig.tight_layout()
    fig.savefig(os.path.join(path_results, 'lollipops', f'{x}.png'), dpi=300)


##
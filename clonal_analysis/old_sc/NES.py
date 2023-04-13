"""
Vizualization of Normalized Enrichment Scores from dist features results
"""

# Code
import sys
import os
import pickle
from Cellula.plotting._plotting_base import *
matplotlib.use('MacOSX')

# Paths
path_main = '/Users/IEO5505/Desktop/MDA_chemo_repro/single_cell/'
path_data = path_main + 'data/'
path_results = path_main + 'results_and_plots/viz/'

# Load data 
nes_1 = pd.read_excel(path_data + 'NES_ENABLE.xlsx', 
    sheet_name='PTs|PC|logit', index_col=0).loc[:, ['NES']]
nes_2 = pd.read_excel(path_data + 'NES_ENABLE.xlsx', 
    sheet_name='lungs|genes|wilcoxon', index_col=0).loc[:, ['NES']]

# The horizontal plot is made using the hline function
nes_1.index = [ ' '.join(x.split(' ')[:-1]) for x in nes_1.index ]
nes_2.index = [ ' '.join(x.split(' ')[:-1]) for x in nes_2.index ]

# Fig
fig, axs = plt.subplots(1,2, figsize=(11,4))
stem_plot(nes_1, 'NES', ax=axs[0])
format_ax(nes_1, ax=axs[0], title='PC8', xlabel='NES', title_size=10, ysize=8)
stem_plot(nes_2, 'NES', ax=axs[1])
format_ax(nes_1, ax=axs[1], title='DE lungs treated vs untreated', xlabel='NES', title_size=10, ysize=8)
fig.tight_layout()
fig.savefig(path_results + 'NES_enable.pdf')
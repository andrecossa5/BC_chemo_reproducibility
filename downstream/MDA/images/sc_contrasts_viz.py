"""
Volcano and lollipop plots of DEGs lists from cusom contrasts.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
from Cellula.plotting._plotting import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'contrasts')


##


# Volcano
df = pd.read_csv(os.path.join(path_results, 'volcano_clones_AACR.csv'), index_colz=0)
fig = volcano(
    df.query('comparison=="g0_vs_g1"'), effect_size='effect_size', evidence='evidence',
    n=5, annotate=False, xlim=(-2, 2)
)
fig.tight_layout()


##
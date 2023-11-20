"""
Final markers.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'general_umaps')

# Read cells meta and format
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata, path_data, contrasts_name='final_cell_state')

# Here we go
D = Dist_features(adata, contrasts, jobs=jobs, app=False)
D.run_all_jobs()

# Save
degs = D.Results.results['cell_state|genes|wilcoxon']['df']
degs.to_csv(os.path.join(path_data, 'DEGs.csv'))


##
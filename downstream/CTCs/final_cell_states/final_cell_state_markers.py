"""
Final markers for final_cell_state.
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *


##


# Paths
path_main = '/Users/ieo7295/Desktop/tests/cell'
path_data = os.path.join(path_main, 'data', 'default')
path_results = "/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/CTCs/general_umaps"

# Read cells meta and format
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))

# Prep contrast and jobs
jobs, contrasts = prep_jobs_contrasts(adata, path_data, contrasts_name='final_cell_state')

# Here we go
D = Dist_features(adata, contrasts, jobs=jobs, app=False)
D.run_all_jobs()

# Save
categories = ['cell_cycle', 'pro-metastatic', 'TNF-alpha_Signaling_via_NF-kB', 
              'Integrated_stress_response', 'interferon','Ubiquitination',
              'mTORC1/folate','oxphos', 'EMT', 'peptidase', 
              'ECM', 'apoptosis', 'DNA_metabolism']
dfs = []

for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = D.Results.results[key]['df'].copy()
    df['label'] = cat  # Add label for source category
    dfs.append(df)

all_degs = pd.concat(dfs)
#degs = D.Results.results['final_cell_state|genes|wilcoxon']['df']
#degs['comparison'].unique()
all_degs.to_csv(os.path.join(path_data, 'DEGs.csv'))


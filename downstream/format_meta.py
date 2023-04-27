#!/usr/bin/python

import sys
import os
import numpy as np
import pandas as pd

# Path meta
path_meta = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/meta/'

# Read cells meta
df_annot = pd.read_csv(os.path.join(path_meta, 'scRNA_batches_200423.csv'), index_col=0)
df = pd.read_csv(os.path.join(path_meta, 'cells_meta.csv'), index_col=0)

# Format
df = (

    df
    .loc[:, ~df.columns.str.contains('passing')]
    .reset_index()
    .merge( 
        df_annot.reset_index().rename(columns={'index':'sample'}),
        on='sample',
        how='left'
    )
    .drop_duplicates()
    .set_index('index')
    .assign(origin=lambda x: np.where(x['sample'].str.contains('PTs'), 'PT', 'lung'))
    .assign(dataset=lambda x: x['sample'].map(lambda x: '_'.join(x.split('_')[:2] + x.split('_')[-1:])) )
    .loc[:, 
        ['GBC', 'sample', 'nUMIs', 'mito_perc', 'detected_genes',\
         'cell_complexity', 'seq_run', 'infection', 'condition', 'origin', 'dataset']
    ]

)

# Save 
df.to_csv(os.path.join(path_meta, 'cells_meta.csv'))
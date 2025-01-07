"""
FIx cells metadata
"""

import os
import pandas as pd
import numpy as np


##


# Args
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'CTCs')

# Save final metadata
meta =  pd.read_csv(os.path.join(path_data, 'cells_meta_orig.csv'), index_col=0)
df = ( 
    meta.groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .to_frame('Prevalence').reset_index()
    .pivot_table(index='GBC', columns='sample', fill_value=0)
)
df.to_csv(os.path.join(path_data, 'clones_albi.csv'))


##


# Fix meta
tests = [meta['sample'].str.contains('PT'), meta['sample'].str.contains('lung'), meta['sample'].str.contains('CTC')] 
meta['origin'] = np.select(tests, ['PT', 'lung', 'CTC'])
meta['dataset'] = meta['sample'].map(lambda x: x.split('_')[1])

# Save meta
meta[['GBC', 'sample', 'dataset', 'origin']+meta.columns[2:-2].to_list()].to_csv(os.path.join(path_data, 'cells_meta.csv'))


##
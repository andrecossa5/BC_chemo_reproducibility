"""
Final cell_assignment MDA dataset.
"""

import os
import pandas as pd
import numpy as np


##


# Args
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'CTCs')

# Save final metadata
meta =  pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)
df = ( 
    meta.groupby('sample')
    ['GBC'].value_counts(normalize=True)
    .to_frame('Prevalence').reset_index()
    .pivot_table(index='GBC', columns='sample', fill_value=0)
)
# df.to_csv(os.path.join(path_data, 'clones_albi.csv'))


##
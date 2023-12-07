"""
Clone calling MDA dataset.
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc


##


# Utils
def summary_what_we_have_longitudinal(meta, branch, dataset, print_top=True):

    PT_sample = f"{branch}_PTs_{dataset}"
    lung_sample = f"{branch}_mets_{dataset}"

    valid_CB_PT = meta.query('sample==@PT_sample').index.map(lambda x: x.split('_')[0])
    valid_CB_lung = meta.query('sample==@lung_sample').index.map(lambda x: x.split('_')[0])

    PT_cells = pd.read_csv(os.path.join(path_clones, PT_sample, 'cells_summary_table.csv'), index_col=0)
    lung_cells = pd.read_csv(os.path.join(path_clones, lung_sample, 'cells_summary_table.csv'), index_col=0)

    valid_CB_PT = set(valid_CB_PT) & set(PT_cells.index)
    valid_CB_lung = set(valid_CB_lung) & set(lung_cells.index)
    valid_clones_PT = PT_cells.loc[list(valid_CB_PT)]['GBC_set'].unique()
    valid_clones_lung = lung_cells.loc[list(valid_CB_lung)]['GBC_set'].unique()

    PT_clones = pd.read_csv(os.path.join(path_clones, PT_sample, 'clones_summary_table.csv'), index_col=0)
    lung_clones = pd.read_csv(os.path.join(path_clones, lung_sample, 'clones_summary_table.csv'), index_col=0)

    set_PT = set(PT_clones.loc[PT_clones['GBC_set'].isin(valid_clones_PT), 'GBC_set'].unique())
    set_lung = set(lung_clones.loc[lung_clones['GBC_set'].isin(valid_clones_lung), 'GBC_set'].unique())

    print(
        f'''
        Dataset {branch}{dataset}: n cells PT {len(valid_CB_PT)}, n cells lung {len(valid_CB_lung)}
        n PT {len(set_PT)}, n lung {len(set_lung)}, common {len(set_PT & set_lung)}
        '''
    )

    PT_clones = PT_clones.set_index('GBC_set').loc[list(set_PT & set_lung)]
    PT_clones.columns = [f'PT {x}' for x in PT_clones.columns]
    lung_clones = lung_clones.set_index('GBC_set').loc[list(set_PT & set_lung)]
    lung_clones.columns = [f'lung {x}' for x in lung_clones.columns]
    top = lung_clones.sort_values('lung prevalence', ascending=False).join(PT_clones)

    if print_top:
        print(top)
        print(top.loc[lambda x: (x['PT n cells']>=10) & (x['lung n cells']>=10)])


##


def what_we_have_all_samples(meta):

    samples = meta['sample'].unique()

    valid_cells = []
    n_cells = {}
    set_clones = {}
   
    for sample in samples:
        valid_CB = meta.query('sample==@sample').index
        cells = pd.read_csv(
            os.path.join(path_clones, sample, 'cells_summary_table.csv'),
            index_col=0
        )
        cells.index = cells.index.map(lambda x: f'{x}_{sample}')
        valid_CB = list(set(valid_CB) & set(cells.index))
        cells = cells.loc[valid_CB].rename(columns={'GBC_set':'GBC'}).join(meta.loc[valid_CB])
        set_clones[sample] = set(cells['GBC'].unique().tolist())
        n_cells[sample] = len(valid_CB)
        valid_cells.append(cells)

    n = len(samples)
    C = np.zeros((n,n))
    for i,x in enumerate(samples):
        for j,y in enumerate(samples):
            C[i,j] = len(set_clones[x] & set_clones[y])

    return (
        pd.DataFrame(C, index=samples, columns=samples), 
        pd.Series(n_cells), 
        pd.concat(valid_cells)
    )


##


# Args
path_data = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/MDA'
path_clones = os.path.join(path_data, 'clonal_info')


##


# Reformat meta: new metadata
meta = pd.read_csv(os.path.join(path_data, 'cells_meta.csv'), index_col=0)

# Add dataset
meta['dataset'] = (
    meta['sample']
    .map(lambda x: '_'.join(x.split('_')[:2]) + '_' + x.split('_')[-1] )
)
meta['origin'] = np.where(meta['sample'].str.contains('mets'), 'lung', 'PT')
L = [
    (meta['sample'].str.startswith('AC_AC')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('AC_AC')) & (meta['origin'] == 'lung'),
    (meta['sample'].str.startswith('NT_NT')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('NT_NT')) & (meta['origin'] == 'lung'),
    (meta['sample'].str.startswith('NT_AC')) & (meta['origin'] == 'PT'),
    (meta['sample'].str.startswith('NT_AC')) & (meta['origin'] == 'lung')
]
c = [
    'PT, treated', 'lung, double-treated', 'PT, untreated', 'lung, untreated',
    'PT, untreated', 'lung, single-treated'
]
meta['condition'] = np.select(L, c)
meta['condition'] = pd.Categorical(
    meta['condition'], 
    categories=[
        'PT, treated', 'PT, untreated', 'lung, double-treated', 
        'lung, single-treated', 'lung, untreated',
    ]
)


##


# Checks single couples, clones
summary_what_we_have_longitudinal(meta, 'AC_NT', 3, print_top=False)


##


# NT_NT: 4+9+3= 16 longitudinal >= 10 cells
# AC_AC: 1+4+x= 5 longitudinal >= 10 cells
# NT_AC: 2+8+x= 5 longitudinal >= 10 cells
# AC_NT: x+x+3= 5 longitudinal >= 10 cells


# Get valid CBC_GBC_sets
common, n_cells, valid_cells = what_we_have_all_samples(meta)
print(n_cells.sum())
print(common.values.mean())


##


# Save final cells
samples = n_cells.loc[lambda x: x>1000].index
valid_cells.to_csv(os.path.join(path_data, 'cells_meta.csv'))

# Format QC.h5ad
adata = sc.read(os.path.join(path_data, 'QC.h5ad'))
adata = adata[meta.index,:].copy()
adata.obs = meta

# Write new adata
adata.write(os.path.join(path_data, 'QC.h5ad'))


##

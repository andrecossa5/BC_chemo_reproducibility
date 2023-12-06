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
    valid_clones_PT = PT_cells.loc[list(valid_CB_PT)]['clone'].unique()
    valid_clones_lung = lung_cells.loc[list(valid_CB_lung)]['clone'].unique()

    PT_clones = pd.read_csv(os.path.join(path_clones, PT_sample, 'clones_summary_table.csv'), index_col=0)
    lung_clones = pd.read_csv(os.path.join(path_clones, lung_sample, 'clones_summary_table.csv'), index_col=0)
    PT_clones = PT_clones.loc[list(valid_clones_PT)]
    lung_clones = lung_clones.loc[list(valid_clones_lung)]

    set_PT = set(PT_clones['GBC_set'].values)
    set_lung = set(lung_clones['GBC_set'].values)

    print(
        f'''
        Dataset {branch}{dataset}: n cells PT {len(valid_CB_PT)}, n cells lung {len(valid_CB_lung)}
        n PT {len(set_PT)}, n lung {len(set_lung)}, common {len(set_PT & set_lung)}
        '''
    )
    print([ ';'.join(x) for x in (set_PT & set_lung)])

    if print_top:
        print('PT')
        print(PT_clones.sort_values('n cells', ascending=False).head(10))
        print('lung')
        print(lung_clones.sort_values('n cells', ascending=False).head(10))


##


def what_we_have_all_samples(meta):

    samples = meta['sample'].unique()

    valid_CBC_GBC_set = []
    n_cells = {}
    set_clones = {}
   
    for sample in samples:

        valid_CB = meta.query('sample==@sample').index.map(lambda x: x.split('_')[0])
        cells = pd.read_csv(os.path.join(path_clones, sample, 'cells_summary_table.csv'), index_col=0)
        valid_CB = set(valid_CB) & set(cells.index)
        valid_clones = cells.loc[list(valid_CB)]['clone'].unique()
        clones = pd.read_csv(os.path.join(path_clones, sample, 'clones_summary_table.csv'), index_col=0)
        clones = clones.loc[list(valid_clones)]
        set_clones[sample] = set(clones['GBC_set'].values)
        n_cells[sample] = len(valid_CB)
        cells_new = (
            cells.loc[list(valid_CB)].reset_index()
            .merge(clones.reset_index(), on='clone')
            [['CBC', 'GBC_set']]
            .set_index('CBC')
        )
        cells_new.index = cells_new.index.map(lambda x: f'{x}_{sample}')
        valid_CBC_GBC_set.append(cells_new)

    n = len(samples)
    C = np.zeros((n,n))
    for i,x in enumerate(samples):
        for j,y in enumerate(samples):
            C[i,j] = len(set_clones[x] & set_clones[y])

    return (
        pd.DataFrame(C, index=samples, columns=samples), 
        pd.Series(n_cells), 
        pd.concat(valid_CBC_GBC_set)
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



##


# Checks single couples, clones
# summary_what_we_have_longitudinal(meta, 'NT_NT', 3, print_top=True)


##


# Get valid CBC_GBC_sets
common, n_cells, cbc_gbc_set = what_we_have_all_samples(meta)
n_cells.sum()
common.values.mean()

# Rename valid GBC_sets
d_names = (
    cbc_gbc_set['GBC_set']
    .value_counts()
    .to_frame('n_cells')
    .assign(name=lambda x: [ f'C{i}' for i in range(x.shape[0]) ])
    ['name'].to_dict()
)
cbc_gbc_set['clone'] = cbc_gbc_set['GBC_set'].map(d_names)
cbc_gbc_set['GBC_set'] = cbc_gbc_set['GBC_set'].astype('str')

# Filter meta with filtered and clonally annotated cells
meta = cbc_gbc_set.join(meta, how='left')


##

# Save final cells
meta.to_csv(os.path.join(path_data, 'cells_meta.csv'))

# Format QC.h5ad
adata = sc.read(os.path.join(path_data, 'QC.h5ad'))
adata = adata[meta.index,:].copy()
adata.obs = meta

# Write new adata
adata.write(os.path.join(path_data, 'QC.h5ad'))


##

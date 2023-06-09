"""
Script to produce pro-metastatic clones summary and visualizations
"""

import os
import numpy as np
import pandas as pd
from Cellula._utils import make_folder
from BC_chemo_utils.plotting import *
from BC_chemo_utils.tests import *


#


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_results = os.path.join(path_main, 'results', 'clonal')

# Read cells meta
df = pd.read_csv(os.path.join(path_main, 'meta', 'cells_meta.csv'), index_col=0)

##

# Here we go
for infection in df['infection'].unique():

    # Filter infection
    df_ = df.query('infection == @infection')
    
    # GBC_sample aggregate, per infection
    df_gbc_sample = (
        df_.groupby(['GBC', 'sample'])
        .size()
        .reset_index(name='count')
        .assign(prevalence=lambda x: \
            x['count'] / x.groupby('sample')['count'].transform('sum')
        )
        .pivot(index='GBC', columns='sample')
        .fillna(0)
    )
    order = ( 
        df_gbc_sample
        .loc[:, pd.IndexSlice['count',:]]
        .sum(axis=1)
        .sort_values(ascending=False)
        .index
    )
    df_gbc_sample.loc[order, :].to_excel(
        os.path.join(path_results, f'df_gbc_sample_infection_{infection}.xlsx')
    )

    test_clone_aggregates(df_, df_gbc_sample)

    # GBC_sample occurrences, per infection
    df_common = (
        df_gbc_sample.loc[:, pd.IndexSlice['count',:]]
        .droplevel(0, axis=1)
        .assign(occurrence=lambda x: np.sum(x>0, axis=1))
        .assign(is_common=lambda x: x['occurrence']>1)
        .sort_values('occurrence', ascending=False)
    )
    df_common.to_excel(
        os.path.join(path_results, f'df_common_infection_{infection}.xlsx')
    )
    
    test_clone_common(df_, df_common)

    # Report
    df_cells_by_sample = df_['sample'].value_counts().to_frame()
    df_infection_clones = df_['GBC'].value_counts().to_frame()
    n_clones = df_infection_clones.shape[0]
    n_clones_10_cells = np.sum(df_infection_clones>10)[0]
    n_common = df_common['is_common'].sum()

    with open(os.path.join(path_results, f'report_{infection}.txt'), 'w') as f:
        f.write(f'# Infection {infection} numbers:\n')
        f.write(f'\n')
        f.write(f'''Total cells and samples: {df_.shape[0]}, 
                    {df_cells_by_sample.shape[0]}\n''')
        f.write(f'n clones: {n_clones}\n')
        f.write(f'n clones > 10 cells {n_clones_10_cells}\n')
        f.write(f'n common: {n_common}\n')
        f.write(f'\n')

    # Save
    df_cells_by_sample.to_excel(
        os.path.join(path_results, f'df_cells_by_sample_{infection}.xlsx'))
    df_infection_clones.to_excel(
        os.path.join(path_results, f'df_clones_total_{infection}.xlsx'))

    # Promet, for each PT, lung couple
    for dataset in df_['dataset'].unique():

        # Folder 
        path_dataset = os.path.join(path_results, dataset)
        make_folder(path_results, dataset)

        # Get samples and subset df_gbc_sample
        samples_names = df_.query('dataset == @dataset')['sample'].unique()
        n_PT = df_.query('dataset == @dataset and origin == "PT"').shape[0]
        n_lung = df_.query('dataset == @dataset and origin == "lung"').shape[0]

        df_promet = (
            df_gbc_sample
            .loc[:, pd.IndexSlice['count', samples_names]]
            .droplevel(0, axis=1)
        )
        df_promet.columns = np.where(df_promet.columns.str.contains('PT'), 'PT', 'lung')
        df_promet = (
            df_promet
            .assign(total=lambda x: x['PT'] + x['lung'])
            .query('PT>0 and lung>0')
            .sort_values(by='total', ascending=False)
            .assign(PT_prevalence=lambda x: x['PT'] / n_PT)
            .assign(lung_prevalence=lambda x: x['lung'] / n_lung)
            .assign(metastatic_potential=lambda x: x['lung_prevalence'] / x['PT_prevalence'])
            .sort_values(by='lung_prevalence', ascending=False)
        )
        
        test_promet(df_, df_promet, dataset)

        # Report
        df_promet.to_excel(os.path.join(path_dataset, f'df_promet.xlsx'))
        with open(os.path.join(path_results, f'report_{infection}.txt'), 'a') as f:
            f.write(f'# Dataset {dataset} pro-metastatic clones: {df_promet.shape[0]}\n')

        # Plot promets and PT clones
        fig, axs = plt.subplots(1,2,figsize=(11,5))
        packed_circle_plot(df_promet, covariate='PT_prevalence', ax=axs[0], 
                           color='orange', annotate=True)
        axs[0].set(title='PT prevalence')
        packed_circle_plot(df_promet, covariate='lung_prevalence', ax=axs[1], 
                           color='r', annotate=True)
        axs[1].set(title='lung prevalence')
        fig.suptitle(f'Dataset {dataset}: metastatic clones')
        fig.savefig(os.path.join(path_dataset, 'promet_circles.png'))


##


# At last write separate excel sheets per sample
writer = pd.ExcelWriter(
    os.path.join(path_results, f'df_gbc_single_samples.xlsx'),
    engine='xlsxwriter'
)

for sample in df['sample'].unique():
    df_ = df.query('sample == @sample')
    gbc_by_sample = (
        df_.groupby('GBC')
        .size()
        .to_frame(name='n_cells')
        .assign(prevalence=lambda x: x['n_cells']/x['n_cells'].sum())
        .sort_values('prevalence', ascending=False)
    )
    gbc_by_sample.to_excel(writer, sheet_name=sample)

writer.save()
"""
Script to produce pro-metastatic clones summary and visualizations
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import circlify


##


def packed_circle_plot(df, covariate=None, ax=None, color='b', annotate=False):
    """
    Circle plot. Packed.
    """

    circles = circlify.circlify(
        df[covariate].tolist(), 
        show_enclosure=True, 
        target_enclosure=circlify.Circle(x=0, y=0, r=1)
    )

    lim = max(
        max(
            abs(c.x) + c.r,
            abs(c.y) + c.r,
        )
        for c in circles
    )
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)

    for name, circle in zip(df.index, circles):
        x, y, r = circle
        ax.add_patch(
            plt.Circle((x, y), r*0.95, alpha=0.4, linewidth=1, 
                fill=True, edgecolor="black", facecolor=color)
        )
        if annotate:
            ax.annotate(name, (x,y), va='center', ha='center', fontsize=4)
    ax.axis('off')

    return ax


##


# Paths
path_main = sys.argv[1]
path_data = os.path.join(path_main, 'data')
path_results = os.path.join(path_main, 'results/clonal_sc')


##


# Read and format QC matrix all samples, meta already formatted
adata = sc.read(os.path.join(path_data, 'QC.h5ad'))
df = adata.obs
df['GBC'] = df['GBC'].astype('str') # Avoid pd.Categorical troubles...
df['sample'] = df['sample'].astype('str')


##


# Here we go
for infection in adata.obs['infection'].unique():

    # GBC_sample aggregate, per infection
    df_ = df.query('infection == @infection')
    df_gbc_sample = (
        df_.groupby(['GBC', 'sample'])
        .size()
        .reset_index(name='count')
        .assign(
            prevalence=lambda x: x
            .groupby('sample')
            .apply(lambda x: x['count'] / x['count'].sum())
            .reset_index(drop=True)
        )
        .pivot(index='GBC', columns='sample')
    )
    df_gbc_sample.to_excel(
        os.path.join(path_results, f'df_gbc_sample_infection_{infection}.xlsx')
    )

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

    # Report
    df_cells_by_sample = df_['sample'].value_counts().to_frame()
    df_infection_clones = df_['GBC'].value_counts().to_frame()
    n_clones = df_infection_clones.shape[0]
    n_clones_10_cells = np.sum(df_infection_clones>10)
    n_common = df_common['is_common'].sum()

    with open(os.path.join(path_results, f'report_{infection}.txt'), 'w') as f:
        f.write(f'# Infection {infection} numbers:\n')
        f.write(f'\n')
        f.write(f'Total cells and samples: {df_.shape[0]}, {df_cells_by_sample.shape[0]}\n')
        f.write(f'n clones: {n_clones}\n')
        f.write(f'n clones > 10 cells {n_clones_10_cells}\n')
        f.write(f'n common: {n_common}\n')
        f.write(f'\n')

    # Save
    df_cells_by_sample.to_excel(os.path.join(path_results, f'df_cells_by_sample_{infection}.xlsx'))
    df_infection_clones.to_excel(os.path.join(path_results, f'df_clones_total_{infection}.xlsx'))

    # Promet, for each PT, lung couple
    for dataset in df['dataset'].unique():

        # Folder 
        path_dataset = os.path.join(path_results, dataset)
        os.mkdir(path_dataset)

        # Get samples and subset df_gbc_sample
        samples_names = df.query('dataset == @dataset')['sample'].unique()

        df_promet = df_gbc_sample.loc[:, 
                    pd.IndexSlice['count', samples_names]].droplevel(0, axis=1)
        df_promet.columns = np.where(df_promet.columns.str.contains('PT'), 'PT', 'lung')
        df_promet = (
            df_promet
            .assign(total=lambda x: x['PT']+x['lung'])
            .query('PT>0 and lung>0')
            .sort_values(by='total', ascending=False)
            .assign(PT_prevalence=lambda x: x['PT']/x['PT'].sum())
            .assign(lung_prevalence=lambda x: x['lung']/x['lung'].sum())
            .assign(metastatic_potential=lambda x: x['lung_prevalence']/x['PT_prevalence'])
            .sort_values(by='metastatic_potential', ascending=False)
        )
        df_promet.to_excel(os.path.join(path_dataset, f'df_promet.xlsx'))

        # Report
        with open(os.path.join(path_results, f'report_{infection}.txt'), 'w') as f:
            f.write(f'# Dataset {dataset} pro-metastatic clones: {df_promet.shape[0]}\n')

        # Plot promets and PT clones
        fig, axs = plt.subplots(1,2,figsize=(10,5))
        packed_circle_plot(df_, covariate='PT_prevalence', ax=axs[0], color='orange', annotate=True)
        axs[0].set(title='PT prevalence')
        packed_circle_plot(df_, covariate='lung_prevalence', ax=axs[1], color='r', annotate=True)
        axs[1].set(title='lung prevalence')
        fig.suptitle(f'Dataset {dataset}: metastatic clones')
        fig.savefig(path_dataset, 'promet_circles.png')

        fig, ax = plt.subplots(figsize=(5,5))
        packed_circle_plot(
            df_gbc_sample.assign(PT_prevalence=lambda x: x['PT']/x['PT'].sum()), 
            covariate='PT_prevalence', ax=ax, color='orange'
        )
        ax.set(title=f'Dataset {dataset}: PT clones')
        fig.savefig(path_dataset, 'PT_circles.png')







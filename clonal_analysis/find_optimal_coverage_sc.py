#!/usr/bin/python

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


##


def cell_assignment(df, t):
    """
    Assing cells to GBCs.
    """
       
    test = (df['read_counts'] > 30) & (df['umi_counts'] > 3) & (df['coverage'] > t)
    df['status'] = np.where(test, 1, 0)

    # Cell MOI tally 
    cell_MOI = (
        df.loc[:, ['CBC', 'GBC', 'status']]
        .query('status == 1')
        .groupby('CBC')
        .sum()
        .rename(columns={'status': 'multiplicity'})
    )
    # Filter uniquely assigned cells and get their unique GBC
    uniquely_assigned_cells = cell_MOI.query('multiplicity == 1').index

    return df_combos, uniquely_assigned_cells


##


# Path
path_sc_results = sys.argv[1]
tresholds = [5, 10, 15, 20, 30, 40, 50, 75]

# Here we go
for sample in os.listdir(path_sc_results):

    if not sample.startswith('.DS_'):

        os.chdir(os.path.join(path_sc_results, sample))
        df_combos = pd.read_csv(
            os.path.join(path_sc_results, sample, 'CBC_GBC_combos.tsv.gz'), 
            index_col=0, 
            sep='\t'
        )

        L = []
        for t in tresholds:
            _, assigned_cells = cell_assignment(df_combos, t)
            L.append(len(assigned_cells))
        optimal_t = tresholds[np.argmax(L)]
        df_combos, uniquely_assigned_cells = cell_assignment(df_combos, optimal_t)

        # Cell assignment plot
        fig, ax = plt.subplots()
        x = np.log10(df_combos['read_counts'])
        y = np.log10(df_combos['umi_counts'])
        ax.plot(x[df_combos['status'] == 1], y[df_combos['status'] == 1], 
                '.', label='assigned', color='blue')
        ax.plot(x[df_combos['status'] == 0], y[df_combos['status'] == 0], 
                '.', label='not-assigned', color='grey')
        ax.set(
            title='CBC-GBC combination status', 
            xlabel='log10_read_counts', 
            ylabel='log10_umi_counts'
        )
        ax.legend()
        fig.savefig(os.path.join(path_sc_results, sample, 'CBC_GBC_combo_status_refined.png'))

        # Compute summary tables: cells 
        df_cells = (
            df_combos
            .query('CBC in @uniquely_assigned_cells and status == 1')
            .drop(columns='status')
            .set_index('CBC')
        )
        
        # Assert we have taken the correct ones and save
        assert df_cells.index.size == uniquely_assigned_cells.size
        df_cells.to_csv(os.path.join(path_sc_results, sample, 'cells_summary_table_refined.csv'))

        # Compute summary tables: clones
        df_clones = df_cells.reset_index().groupby('GBC').size().to_frame('n_cells')
        df_clones['prevalence'] = df_clones['n_cells'] / df_clones['n_cells'].sum()
        df_clones = df_clones.sort_values('prevalence', ascending=False)
        df_clones.to_csv(os.path.join(path_sc_results, sample, 'clones_summary_table_refined.csv'))
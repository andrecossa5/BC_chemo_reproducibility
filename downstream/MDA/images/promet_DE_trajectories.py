"""
Patterns of promet contrasts (DE) and clonal trajectories
"""

import os
import numpy as np
import pandas as pd
from plotting_utils._plotting_base import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'MDA', 'DE_traj')

# Read gene lists
df_nt = pd.read_csv(os.path.join(path_data, 'promet_NT_DEGs.csv'), index_col=0)
df_nt = df_nt.loc[df_nt['comparison']=='g0_vs_g1']
df_ac = pd.read_csv(os.path.join(path_data, 'promet_AC_DEGs.csv'), index_col=0)
df_ac = df_ac.loc[df_ac['comparison']=='g0_vs_g1']
nt_traj = pd.read_csv(os.path.join(path_data, 'NT_NT_met_potential_traj.csv'), index_col=0)
ac_traj = pd.read_csv(os.path.join(path_data, 'AC_AC_met_potential_traj.csv'), index_col=0)

# Merge
df_de = (
    df_nt.rename(columns={'evidence':'FDR', 'effect_size':'log2FC'})[['FDR', 'log2FC','AUROC']]
    .join(
        df_ac.rename(columns={'evidence':'FDR', 'effect_size':'log2FC'})[['FDR', 'log2FC','AUROC']],
        lsuffix='_NT', rsuffix='_AC'
    )
    .assign(mean_AUROC=lambda x: x[['AUROC_NT', 'AUROC_AC']].mean(axis=1))
)
df_traj = nt_traj.join(ac_traj, lsuffix='_NT', rsuffix='_AC')


##


# Save and plot relevant
fig, ax = plt.subplots(figsize=(4,4))
sns.kdeplot(df_de['mean_AUROC'], ax=ax, fill=True, alpha=.3)
ax.axvline(x=.6, color='r', linestyle='--', linewidth=1)
n_selected = np.sum(df_de['mean_AUROC']>=.6)
ax.text(.65,.9,f'n genes: {n_selected}',transform=ax.transAxes)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'distr_auroc.png'), dpi=300)


# Save
# pd.Series(df_de.query('mean_AUROC>=.6').index).to_csv(os.path.join(path_data, 'top_genes_promet_DEGs_AUROC.csv'))


##


# Viz importances: DE
fig, axs = plt.subplots(1,3,figsize=(9,3))
axs[0].plot(df_de['log2FC_NT'],df_de['log2FC_AC'],'ko',markersize=1)
axs[0].set(xlabel='log2FC, NT', ylabel='log2FC AC')
axs[1].plot(-np.log10(df_de['FDR_NT']),-np.log10(df_de['FDR_AC']),'ko',markersize=1)
axs[1].set(xlabel='-log10(FDR), NT', ylabel='-log10(FDR), AC')
axs[2].plot(df_de['AUROC_NT'],df_de['AUROC_AC'],'ko',markersize=1)
axs[2].set(xlabel='AUROC, NT', ylabel='AUROC, AC')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'DE_contrasts.png'), dpi=300)

# Viz AUROC
fig, ax = plt.subplots(figsize=(4,4))
ax.plot(df_de['AUROC_NT'],df_de['AUROC_AC'],'ko',markersize=1)
ax.set(xlabel='AUROC, NT', ylabel='AUROC, AC')
ax.axvline(x=.5, color='red', linestyle='--', linewidth=1)
ax.axhline(y=.5, color='red', linestyle='--', linewidth=1)
test = (df_de['AUROC_NT']>.63) & (df_de['AUROC_AC']>.63)
ta.allocate_text(
    fig, ax, 
    df_de.loc[test]['AUROC_NT'], df_de.loc[test]['AUROC_AC'], df_de.loc[test].index,
    x_scatter=df_de['AUROC_NT'], y_scatter=df_de['AUROC_AC'], 
    linecolor='k', textcolor='r', textsize=8, 
    max_distance=0.1, linewidth=0.5, nbr_candidates=100
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'DE_AUROC.png'), dpi=300)


##


# Viz importances: trajectories
fig, axs = plt.subplots(1,2,figsize=(6,3))
axs[0].plot(df_traj['slopes_NT'],df_traj['slopes_AC'],'ko',markersize=1)
axs[0].set(xlabel='Slopes, NT', ylabel='Slopes, AC')
axs[1].plot(-np.log10(df_traj['qs_NT']),-np.log10(df_traj['qs_AC']),'ko',markersize=1)
axs[1].set(xlabel='-log10(qvalue), NT', ylabel='-log10(qvalue), AC')
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'traj.png'), dpi=300)

# Viz slopes
fig, ax = plt.subplots(figsize=(4,4))
ax.plot(df_traj['slopes_NT'],df_traj['slopes_AC'],'ko',markersize=1)
ax.set(xlabel='Slopes, NT', ylabel='Slopes, AC')
ax.axvline(x=0, color='red', linestyle='--', linewidth=1)
ax.axhline(y=0, color='red', linestyle='--', linewidth=1)
test = (df_traj['slopes_NT']>.05) & (df_traj['slopes_NT']>.05)
ta.allocate_text(
    fig, ax, 
    df_traj.loc[test]['slopes_NT'], df_traj.loc[test]['slopes_AC'], df_traj.loc[test].index,
    x_scatter=df_traj['slopes_NT'], y_scatter=df_traj['slopes_AC'], 
    linecolor='k', textcolor='r', textsize=8, 
    max_distance=0.1, linewidth=0.5, nbr_candidates=100
)
fig.tight_layout()
fig.savefig(os.path.join(path_results, 'traj_slopes.png'), dpi=300)


##



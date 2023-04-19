"""
Vizualization of categorical covariates plus single_cell clones in separate datasets
"""

# Code
import os
import sys
import Cellula
from Cellula.preprocessing._embeddings import * 
from Cellula.plotting._plotting import *
from Cellula._utils import make_folder
matplotlib.use('MacOSX')

# Paths
path_main = '/Users/IEO5505/Desktop/MDA_chemo_repro/single_cell/'
make_folder(path_main + 'results_and_plots/viz/', 'clonal')
make_folder(path_main + 'results_and_plots/report/', 'clonal')
path_report = path_main + 'results_and_plots/report/clonal/'
path_viz = path_main + 'results_and_plots/viz/clonal/'


##


# Data
adata = sc.read(path_main + 'data/clustered.h5ad') 
adata.obs['origin'] = adata.obs['origin'].astype(str)
adata.obs['origin'][adata.obs['origin'] == 'nan'] = 'reference'
adata.obs['origin'] = pd.Categorical(adata.obs['origin'])
adata.obs['dataset'] = adata.obs['dataset'].astype(str)
adata.obs['dataset'][adata.obs['dataset'] == '0'] = 'reference'
adata.obs['dataset'] = pd.Categorical(adata.obs['dataset'])

# Filter at least 10 cells in a sample
ncells = 10
GBCs = adata.obs.groupby(['GBC', 'sample']).size() # ... total clones
GBCs_filtered = GBCs[GBCs>ncells].reset_index('GBC')['GBC'].astype(str).unique() # n clones with more than 10 cells

# Subset cells
adata = adata[adata.obs.query('GBC in @GBCs_filtered').index, :] # adata.obs['GBC']: pd.Series with cat dtype
meta = adata.obs

##


################################################################

# Clones counts

# Utils
def cc_counts(meta, cov):
    """
    Count clones and cells per .obs grouping.
    """
    clones_counts = []
    cells_counts = []
    samples = meta[cov].unique()
    for s in samples:
        n_clones = np.sum(meta.groupby([cov, 'GBC']).size()[s] > 0)
        n_cells = meta.query(f'{cov} == @s').shape[0]
        clones_counts.append(n_clones)
        cells_counts.append(n_cells)
    df = pd.DataFrame(
        {'n_cells' : cells_counts, 'n_clones' : clones_counts}, 
        index=samples
    ).sort_values('n_clones', ascending=False)

    return df

##

def plot_counts_one(meta, cov):
    fig, axs = plt.subplots(1, 2)
    df_ = cc_counts(meta, cov)
    bar(df_, y='n_cells', ax=axs[0], c='lightblue', s=0.6)
    format_ax(df_, axs[0], title='n_cells', ylabel='n', xticks=df_.index, rotx=90)
    bar(df_, y='n_clones', ax=axs[1], c='blue', s=0.6)
    format_ax(df_, axs[1], title='n_clones', ylabel='n', xticks=df_.index, rotx=90)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    return fig

# Create and save tables and viz
with PdfPages(path_viz + f'cells_and_clones_numbers.pdf') as pdf:
    for cov in [ 'sample', 'condition', 'dataset', 'cc_phase', 'origin']: 
        cc_counts(adata.obs, cov).to_excel(path_report + f'{cov}.xlsx')
        fig = plot_counts_one(meta, cov)
        fig.suptitle(cov.capitalize())
        pdf.savefig()  
        plt.close()


##


################################################################

# Dataset clones stats

# Clones stats
def dataset_clones_stats(adata, dt, by='FC'):
    """
    Calculate some important clone stats, by dataset.
    """
    meta = adata.obs.astype({'GBC' : 'str'})

    meta['dataset'].unique()

    df_ = meta.query('dataset == @dt')

    n_c_PT = []
    f_PT = []
    n_c_lung = []
    f_lung = []
    n_c_ref = []
    f_ref = []
    is_common = []
    in_ref = []

    x = df_.query('origin == "PT"')['GBC'].unique()
    y = df_.query('origin == "lung"')['GBC'].unique()
    z = meta.query('origin == "reference"')['GBC'].unique()

    for clone in df_['GBC'].unique():
        d_clone = df_.query('GBC == @clone')
        n_cells_PT = d_clone.query('origin == "PT"').shape[0]
        n_cells_lung = d_clone.query('origin == "lung"').shape[0]
        n_cells_ref = meta.query('origin == "reference" & GBC == @clone').shape[0]
        commm = (clone in x) and (clone in y)
        refff = ((clone in x) and (clone in z)) or ((clone in x) and (clone in z))
        
        n_c_PT.append(n_cells_PT)
        f_PT.append(n_cells_PT / df_.query('origin == "PT"').shape[0])
        n_c_lung.append(n_cells_lung)
        f_lung.append(n_cells_lung / df_.query('origin == "lung"').shape[0])
        n_c_ref.append(n_cells_ref)
        f_ref.append(n_cells_ref / meta.query('origin == "reference"').shape[0])
        is_common.append(commm)
        in_ref.append(refff)
    
    df = pd.DataFrame(
        {
            'n_cells_PT' : n_c_PT,
            'f_PT' : f_PT,
            'n_cells_lung' : n_c_lung,
            'f_lung' : f_lung,
            'n_cells_ref' : n_c_ref,
            'f_ref' : f_ref, 
            'is_common' : is_common,
            'in_ref' : in_ref
        }, 
        index=df_['GBC'].unique()
    )

    df['FC'] = df['f_lung'] / (df['f_PT'] + 0.000001)
    df.loc[~df['is_common'], ['FC']] = 0.0
    df = df.sort_values(by=by, ascending=False)

    return df


#


# Create and save tables and viz
datasets = [ x for x in adata.obs['dataset'].cat.categories if x != 'reference' ]
df_ = pd.concat(
    [
        dataset_clones_stats(adata, dt, by='f_lung').query('is_common == True').assign(dataset=dt) \
        for dt in datasets
    ], 
    axis=0
)
df_ = df_.loc[:, 
    [ 'dataset', 'n_cells_ref', 'n_cells_PT', 'n_cells_lung', 'f_PT', 'f_lung', 'FC', 'f_ref', 'is_common', 'in_ref' ]
]
df_.to_excel(path_report + 'pro_metastatic_stats.xlsx')

df_.query('n_cells_ref > 10 and n_cells_PT > 10 and n_cells_lung > 10') 
# No clones with at least 10 cells at each timepoint, within the same dataset
df_.query('n_cells_PT > 10 and n_cells_lung > 10').to_excel(path_report + 'possible_clones_trajectories.xlsx') 
# 12 if only PT and lung. 2 in treated, 10 untreated dataset. If n cells is raised at 50, only 4, one per dataset.

# Viz
with PdfPages(path_viz + 'pro_metastatic_stats.pdf') as pdf:
    for dt in df_['dataset'].unique(): 

        # Data
        df_dt = df_.query('dataset == @dt')

        # Fig
        fig, axs = plt.subplots(1, 3, figsize=(10, 5))

        d = df_dt.loc[:, ['n_cells_PT', 'n_cells_lung']]
        axs[0].plot(df_dt.loc[:, ['n_cells_PT', 'n_cells_lung']].T.values, '--o', c='lightblue')
        axs[0].set_yscale('log')
        axs[0].set_xticks(np.arange(0, len(d.columns), 1))
        axs[0].set_xticklabels(['n_cells_PT', 'n_cells_lung'])
        axs[0].set(title='n_cells', xlabel='', ylabel='value')

        d = df_dt.loc[:, ['f_PT', 'f_lung']]
        axs[1].plot(d.T.values, '--o', c='blue')
        axs[1].set_xticks(np.arange(0, len(d.columns), 1))
        axs[1].set_xticklabels(['f_PT', 'f_lung'])
        axs[1].set(title='frequency', xlabel='', ylabel='')
        axs[1].text(.7, .3, f'n_clones: {df_dt.shape[0]}')

        g = sns.stripplot(
            data=df_dt.loc[:, ['FC', 'f_lung']].rename(columns={'FC':'value'}).assign(metric='FC'), 
            x='metric', y='value', hue='f_lung', palette='viridis', ax=axs[2]
        )
        axs[2].set(title='Fold Change', xlabel='', ylabel='')

        if dt in ['C1_3_sign', 'C1_Null']: 
            g.legend_.remove()

        fig.suptitle(dt.capitalize())
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        pdf.savefig()  
        plt.close()


##


################################################################

# Dataset condition clones stats
meta = adata.obs.astype({'GBC' : 'str', 'condition' : 'str', 'sample' : 'str'})
cond_df = meta.groupby(['GBC', 'condition', 'sample']).size().reset_index(['condition', 'sample']).rename(
    columns={0:'n_cells'}
)
cond_df = cond_df.query('n_cells > 10 and condition != "ref"')

# Treated PT
PT_df = cond_df.loc[cond_df['sample'].str.startswith('PT')].drop_duplicates()
common_clones_treated_PT_PTs = PT_df.query('condition == "treated_PT"').reset_index('GBC').groupby('GBC').size()
common_clones_treated_PT_PTs = common_clones_treated_PT_PTs[common_clones_treated_PT_PTs > 1].index

meta.query('GBC in @common_clones_treated_PT_PTs').groupby(['GBC', 'sample']).size().to_excel(
    path_report + 'common_clones_PT.xlsx'
)
# 3 common clones in the treated PT branch which have more than 10 cells in their PTs

lung_df = cond_df.loc[cond_df['sample'].str.startswith('lung')].drop_duplicates()
common_clones_treated_PT_lungs = lung_df.query('condition == "treated_PT"').reset_index('GBC').groupby('GBC').size()
common_clones_treated_PT_lungs = common_clones_treated_PT_lungs[common_clones_treated_PT_lungs > 1].index
# No common clones in the treated PT branch which have more than 10 cells in their lungs

# Untreated PT
PT_df = cond_df.loc[cond_df['sample'].str.startswith('PT')].drop_duplicates()
common_clones_untreated_PT_PTs = PT_df.query('condition == "untreated_PT"').reset_index('GBC').groupby('GBC').size()
common_clones_untreated_PT_PTs = common_clones_untreated_PT_PTs[common_clones_untreated_PT_PTs > 1].index

meta.query('GBC in @common_clones_untreated_PT_PTs').groupby(['GBC', 'sample']).size() 
# 1 common clone in the untreated PT branch which have more than 10 cells in its PTs

lung_df = cond_df.loc[cond_df['sample'].str.startswith('lung')].drop_duplicates()
common_clones_untreated_PT_lungs = lung_df.query('condition == "untreated_PT"').reset_index('GBC').groupby('GBC').size()
common_clones_untreated_PT_lungs = common_clones_untreated_PT_lungs[common_clones_untreated_PT_lungs > 1].index
meta.query('GBC in @common_clones_untreated_PT_lungs').groupby(['GBC', 'sample']).size() 
# 1 common clones in the untreated PT branch which have more than 10 cells in their lungs

# Across condition
for x in meta['GBC'].unique():
    d_ = meta.query('GBC == @x')
    d_ = d_.loc[d_['sample'].str.startswith('lung')] # lung
    if d_.shape[0] > 0:
        occ_samples = d_.groupby('sample').size()
        if np.sum(occ_samples>0) > 1:
            print(x)
            print(occ_samples)

meta.query('GBC == "TATAGGTGTACGTTACGG"').groupby('sample').size()
# 1 common clone across conditions in PTs (0 in lungs): TATAGGTGTACGTTACGG, altri 'recuperabili... 1 2 forse.

################################################################

"""
Clonal utils.
"""

import numpy as np
import pandas as pd


##


def stats_summary(df, key='sample', freq='freq_sc'):
    """
    Basic stats summary.
    """
    stats = {
        key : [],
        'n' : [],
        'min_freq' : [],
        'max_freq' : [],
        'median_freq' : [],
        'mean_freq' : [],
        'std_freq' : [],
        'SH' : [],
        'EI' : []
    }
    
    for c in df['sample'].unique():

        df_one = df.query('sample == @c')
        stats_one = df_one[freq].describe()

        stats[key].append(c)
        stats['n'].append(df_one.shape[0])
        stats['min_freq'].append(stats_one['min'])
        stats['max_freq'].append(stats_one['max'])
        stats['mean_freq'].append(stats_one['mean'])
        stats['median_freq'].append(stats_one['50%'])
        stats['std_freq'].append(stats_one['std'])
        stats['SH'].append(SH(df_one[freq]))
        stats['EI'].append(EI(df_one[freq]))

    return pd.DataFrame(stats).set_index(key)


##


def SH(p):
    """
    Shannon entropy of some pdist p.
    """
    return -np.sum( np.log10(p) * p )


##


def EI(p):
    """
    Expansion Index of some pdist p.
    """
    return 1 -  ( SH(p) / np.log10(p.size) )


## 


def common(df):
    """
    Median and mean common clones across samples, plus their number for all pairwise comparisons
    """

    samples = df['sample'].unique()
    n_samples = samples.size

    counts = np.zeros((n_samples, n_samples))
    common_idxs = {}

    for i, x in enumerate(samples):
        GBC_x = set(df.query('sample == @x').index.to_list())
        for j, y in enumerate(samples):
            if j>=i:
                GBC_y = set(df.query('sample == @y').index.to_list())
                inter =  GBC_x & GBC_y
                common_idxs[f'{x}|{y}'] = inter
                counts[i, j] = len(inter)

    d = { 'median' : np.median(counts), 'mean' : np.mean(counts) }

    return d, counts, common_idxs
    

##


def common_in(df, n, normalize=False):
    """
    % of clones shared by at least n samples.
    """

    GBCs = df.index.unique()
    common_counts = []

    for x in GBCs:
        d = df.loc[x, 'sample']
        try:
            common_counts.append(df.loc[x, 'sample'].shape[0])
        except:
            common_counts.append(1)

    if normalize:
        return np.sum(np.array(common_counts) >= n) / GBCs.size
    else:
        return np.sum(np.array(common_counts) >= n) 


#


def calculate_clonal_SH_and_prevalence(adata, cov='cell_states'):
    """
    Util to calculate for each clone-sample combo the SH of its cell state frequency and its cell state prevalences.
    """

    df_ = (
        adata.obs.groupby(['condition', 'origin', 'sample'])
        ['GBC'].value_counts(normalize=True).loc[lambda x:x>0]
        .reset_index(name='freq')
        .rename(columns={'level_3':'GBC'})
        .merge(
            adata.obs[['sample', 'dataset']]
            .reset_index(drop=True)
            .drop_duplicates(),

        )
        .drop_duplicates()
        .pivot_table(index='GBC', columns=['dataset', 'origin'], values='freq')
        .melt(value_name='freq', ignore_index=False)
        .reset_index()
        .pivot_table(index=['GBC', 'dataset'], columns='origin', values='freq')
        .reset_index().dropna().set_index('GBC')
        .assign(met_potential=lambda x: x['lung']/x['PT'])
    )

    # SH and frequency of each clone in a cluster
    SH_PT = []
    SH_lung = []
    P_L = []
    for i in range(df_.shape[0]):

        GBC = df_.index[i]
        dataset = df_.iloc[i,0]
        SH_PT.append(
        SH(
            adata.obs
            .query('GBC==@GBC and dataset==@dataset and origin=="PT"')
            [cov].value_counts(normalize=True).values + 0.00000001
        ))
        SH_lung.append(
        SH(
            adata.obs.query('GBC==@GBC and dataset==@dataset and origin=="lung"')
            [cov].value_counts(normalize=True).values + 0.00000001
        ))
        P_L.append(
            adata.obs.query('GBC==@GBC and dataset==@dataset')[cov]
            .value_counts(normalize=True)
            .loc[adata.obs[cov].cat.categories]
            .values[np.newaxis,:]
        )

    df_prev = pd.DataFrame(
        np.concatenate(P_L, axis=0), 
        columns=adata.obs[cov].cat.categories,
        index=df_.index
    )

    df_['SH_PT'] = SH_PT
    df_['SH_lung'] = SH_lung
    df_['diff_SH'] = df_['SH_PT'] / df_['SH_lung']
    df_ = df_.join(df_prev)

    return df_


##
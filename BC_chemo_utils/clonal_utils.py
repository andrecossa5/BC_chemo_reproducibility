"""
Clonal utils.
"""

import numpy as np
import pandas as pd


##


def stats_summary(df, key='sample', use_spikeins=False):
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
    
    for c in df['sample'].cat.categories:

        df_one = df.query('sample == @c')
        freq = 'cellular_prevalence_wi' if use_spikeins else 'cellular_prevalence_wo'
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

    samples = df['sample'].cat.categories
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

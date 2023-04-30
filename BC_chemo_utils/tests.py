"""
Test functions.
"""

import numpy as np
import pandas as pd


#


def test_clone_aggregates(df_, df_gbc_sample):

    df_counts = (
        df_gbc_sample
        .loc[:, pd.IndexSlice['count',:]]
        .droplevel(None, axis=1)
    )
    df_freq = (
        df_gbc_sample
        .loc[:, pd.IndexSlice['prevalence',:]]
        .droplevel(None, axis=1)
    )
    
    assert df_counts.values.sum() == df_.shape[0]
    assert df_freq.values.min() <= 1
    assert df_freq.values.min() >= 0
    
    for s in df_freq.columns:
        clones = df_freq.index[np.random.choice(df_freq.index.size, 10)]
        for clone in clones:
            count_sample = df_.query('sample == @s').shape[0]
            count_clone = df_.query('sample == @s and GBC == @clone').shape[0]
            
            assert count_clone / count_sample == df_freq.loc[clone, s]
            
            
##


def test_clone_common(df_, df_common):

    clones = df_common.index[np.random.choice(df_common.index.size, 10)]
    samples = df_common.columns[:-2]
    
    for c in clones:
        occurrence = df_common.loc[c,'occurrence']
        L = []
        for s in samples:
            L.append(c in df_.query('sample == @s')['GBC'].values)
        
        assert occurrence == np.array(L).sum() 


##


def test_promet(df_, df_promet, dataset):
    
    df_dt = df_.query('dataset == @dataset')
    clones = df_promet.index[:10] if df_promet.shape[0] > 10 else df_promet.index
    
    for c in clones:
        PT,  lung,  met = df_promet.loc[c][-3:]
        
        n_clone_PT = df_dt.query('GBC == @c and origin == "PT"').shape[0]
        n_PT = df_dt.query('origin == "PT"').shape[0]
        n_clone_lung = df_dt.query('GBC == @c and origin == "lung"').shape[0]
        n_lung = df_dt.query('origin == "lung"').shape[0]
        
        assert PT == (n_clone_PT / n_PT)
        assert lung == (n_clone_lung / n_lung)
        assert met == (lung / PT)       
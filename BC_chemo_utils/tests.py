"""
Test functions.
"""

import numpy as np
import pandas as pd
from gseapy import prerank
from joblib import cpu_count
from scipy.special import binom
from scipy.stats import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_percentage_error, r2_score



##


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


##


def compute_enrichment(df, col1, col2, target):
    """
    Compute -log10(FDR) Fisher's exact test: col1 enrichment towards col2 target category.
    """

    n = df.shape[0]
    groups = np.sort(df[col1].unique())

    target_ratio_array = np.zeros(groups.size)
    oddsratio_array = np.zeros(groups.size)
    pvals = np.zeros(groups.size)

    # Here we go
    for i, g in enumerate(groups):

        test_group = df[col1] == g
        test_target = df[col2] == target

        group_size = test_group.sum()
        group_target_size = (test_group & test_target).sum()
        target_ratio = group_target_size / group_size
        target_ratio_array[i] = target_ratio
        other_groups_state_size = (~test_group & test_target).sum()

        # Fisher
        oddsratio, pvalue = fisher_exact(
            [
                [group_target_size, group_size - group_target_size],
                [other_groups_state_size, n - other_groups_state_size],
            ],
            alternative='greater',
        )
        oddsratio_array[i] = oddsratio
        pvals[i] = pvalue

    # Correct pvals --> FDR
    pvals = multipletests(pvals, alpha=0.05, method="fdr_bh")[1]

    # Results
    results = pd.DataFrame({
        'group': groups,
        'perc_in_target' : target_ratio_array,
        'odds_ratio' : oddsratio_array,
        'FDR' : pvals,
        'enrichment' : -np.log10(pvals) 
    }).sort_values('enrichment', ascending=False)

    # Reformat FDR
    results['FDR'] = results['FDR'].map(lambda x: f'{x:.2e}')

    return results


##


def fastGSEA(s, collection='GO_Biological_Process_2021', n_top=50):
    """
    Quick GSEA with gseapy.
    """

    results = prerank(
        rnk=s,
        gene_sets=[collection],
        threads=cpu_count(),
        min_size=15,
        max_size=500,
        permutation_num=200, 
        outdir=None, 
        seed=1234,
        verbose=True,
    )

    filtered_df = (
        results.res2d.loc[:, [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]]
        .rename(columns={'FDR q-val' : 'p_adjusted'})
        .query('p_adjusted<=.1')
        .sort_values('NES', ascending=False)
        .iloc[:n_top, :]
    )
    pd.options.mode.chained_assignment = None # Remove warning
    new_term = filtered_df['Term'].map(lambda x: x.split('__')[1])
    filtered_df.loc[:, 'Term'] = new_term
    filtered_df = filtered_df.set_index('Term')

    return filtered_df


##


def simple_lm(df, y='met_potential', cov=None, n=100, prop_test=.2):

    MAPEs = []
    R2 = []
    w = []
    for _ in range(n):

        model = LinearRegression(fit_intercept=False)
        scaler = StandardScaler()
        train, test = train_test_split(df, test_size=prop_test)

        X_train_scaled = scaler.fit_transform(train[cov])
        y_train_scaled = scaler.fit_transform(train[y].values.reshape(-1,1))
        X_test_scaled = scaler.fit_transform(test[cov])
        y_test_scaled = scaler.fit_transform(test[y].values.reshape(-1,1))

        model.fit(X_train_scaled, y_train_scaled)
        MAPEs.append(
            mean_absolute_percentage_error(
            model.predict(X_test_scaled), scaler.fit_transform(y_test_scaled) 
        ))
        R2.append(
            r2_score(
            model.predict(X_test_scaled), scaler.fit_transform(y_test_scaled) 
        ))
        w.append(
            pd.Series(model.coef_.flatten(), index=cov).sort_values(ascending=False)
        )
        df_res = pd.concat(w, axis=1)
        order = df_res.median(axis=1).sort_values(ascending=False).index
        df_res = df_res.loc[order]

    return MAPEs, R2, df_res


##
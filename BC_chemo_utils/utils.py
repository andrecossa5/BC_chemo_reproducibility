"""
Miscellanea.
"""

import os
import re
import pandas as pd


##


def read_model_results(path_results, model='fixed'):

    L = []
    for path_root, _, files in os.walk(os.path.join(path_results, model)):
        for x in files:
            L.append(os.path.join(path_root, x))
    df_l = []
    for x in L:
        df_l.append(
            pd.read_csv(x, index_col=0)
            .assign(regressed_cc=bool(re.search('G2.M', x)))
            # .sort_values('FDR', ascending=True)
        )
    df = pd.concat(df_l)

    return df


##


def read_GSEA(path_results, collection='GO_biological_process_2021'):
    
    L = []
    for path_root, _, files in os.walk(os.path.join(path_results, 'GSEA', collection)):
        for x in files:
            splitted = x.split('_')
            contrast = '_'.join(splitted[:-2])
            model = splitted[-2]
            cc = splitted[-1].split('.')[0]
            L.append(
                pd.read_csv(os.path.join(path_root, x), index_col=0)
                .assign(contrast=contrast, model=model, cc=cc)
            )
    df = pd.concat(L)

    return df

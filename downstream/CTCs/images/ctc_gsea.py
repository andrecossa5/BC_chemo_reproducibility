"""
script for GSEA plot
"""

import os
import textalloc as ta
from typing import Dict, Iterable, Any, Tuple
from matplotlib.lines import Line2D 
import anndata as ad
import pandas as pd
import scanpy as sc
import seaborn as sns
import pickle
import time
import gc
import re
from gseapy import prerank
from joblib import cpu_count
import matplotlib.pyplot as plt
from Cellula.dist_features._dist_features import *
from Cellula.dist_features._Dist import *
from Cellula.dist_features._Gene_set import *
import plotting_utils as plu

#Paths
path_main = "/Users/ieo7295/Desktop/BC_chemo_reproducibility"
path_data= os.path.join(path_main,"data", "CTCs") #"grn","resources"
path_results= os.path.join(path_main, "results", "CTCs", "gene_reg_net")

#Data
adata= sc.read_h5ad(os.path.join(path_data, "clustered.h5ad"))
with open(os.path.join(path_data,'pro_CTC_contrasts.pickle'),'rb') as f:
    pro_ctc = pickle.load(f)

with open(os.path.join(path_data,'pt_ctc.pickle'), 'rb') as p:
    pt_ctc = pickle.load(p)

with open(os.path.join(path_data,'pro_CTC_all.pickle'), 'rb') as t:
    pro_ctc_all = pickle.load(t)

with open(os.path.join(path_data,'pro_CTC_all_new.pickle'), 'rb') as s:
    pro_ctc_all_new = pickle.load(s)

with open(os.path.join(path_data,'pro_lung_pt.pickle'), 'rb') as s:
    pro_lung_pt = pickle.load(s)
## DE ##

# Prep contrast and jobs
dfs = []
categories= ['1_proCTC_rest','2_proCTC_PT_rest','3_proCTC_PT_rest','4_proCTC_rest']
category = ['pro_CTC_all']
categ = ['pt_ctc','lung_ctc','pt_lung','ctc_pt/lung']
gories = ['pro_PT_vs_CTC','pro_lung_vs_CTC','pro_vs_rest_PT','pro_vs_rest_CTC','pro_vs_rest_lung']
cate = ['pro_lung_vs_PT']


for cat in categories:
    key = f"{cat}|genes|wilcoxon"
    df = pro_ctc.results[key]['df'].copy()
    df['label'] = cat 
    dfs.append(df)

for cat in category:
    key = f"{cat}|genes|wilcoxon"
    df = pro_ctc_all.results[key]['df'].copy()
    df['label'] = cat
    dfs.append(df)

for cat in categ:
    key = f"{cat}|genes|wilcoxon"
    df = pt_ctc.results[key]['df'].copy()
    df['label'] = cat
    dfs.append(df)

for cat in gories:
    key = f"{cat}|genes|wilcoxon"
    df = pro_ctc_all_new.results[key]['df'].copy()
    df['label'] = cat
    dfs.append(df)

for cat in cate:
    key = f"{cat}|genes|wilcoxon"
    df = pro_lung_pt.results[key]['df'].copy()
    df['label'] = cat
    dfs.append(df)

def fastGSEA(s, collection='MSigDB_Hallmark_2020', n_top=50):  #GO_Biological_Process_2023 #MSigDB_Hallmark_2020
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

all_degs= pd.concat(dfs)
all_degs.to_csv(os.path.join(path_data, 'CTC_DEGs_gsea.csv'))
already_done= ['1_proCTC_rest','2_proCTC_PT_rest','pro_CTC_all','3_proCTC_PT_rest','4_proCTC_rest','pt_ctc',  
               'lung_ctc','pt_lung','ctc_pt/lung', 'pro_PT_vs_CTC','pro_lung_vs_CTC',                     
               'pro_vs_rest_PT','pro_vs_rest_CTC','pro_vs_rest_lung']
for i, df in enumerate(dfs):
    if df['label'].iloc[0] == 'pro_lung_vs_PT':
        try:
            rank= df.query('comparison == "g1_vs_g0"')['effect_size'].sort_values(ascending=False)
            name_file = df['label'].iloc[0]
            df_ = fastGSEA(rank)
            df_.to_csv(
                os.path.join(
                    path_results, 'pro_lung_pt',f"{name_file}_g1_vs_g0_Hallmark"
                )
            )
            time.sleep(2)
            del rank, df_, name_file
            gc.collect()
        except Exception as e:
            print(f"Failed on index {i} ({df['label'].iloc[0] if 'label' in df else 'unknown'}): {e}")



# Plot     
for path_root, _, files in os.walk(os.path.join(path_results, 'pro_lung_pt')): #GO_biological_2023
    for x in files:
        if x.endswith('0_Hallmark'):
            df = pd.read_csv(os.path.join(path_root, x), index_col=0)
            fig, ax = plt.subplots(figsize=(8.5,5))
            plu.stem_plot(df.head(20), 'NES', ax=ax)
            plu.format_ax(ax, xlabel='NES', yticks=df.head(20).index.map(lambda x: x.split('(')[0] ), title='Top enriched pathways')
            fig.tight_layout()
            fig.savefig(os.path.join(path_results, 'pro_lung_pt', f'{x}_Hallmark.png'), dpi=300)
            plt.close(fig)


## AA transporter analysis
aa_df = pd.read_csv(os.path.join(path_data,'aa_transporters.txt')).drop(columns=['Unnamed: 1'])
aa_list = aa_df['Transporter'].tolist()
elchain_genes = ['ATP5F1A','UQCRC2','SDHB','NDUFB8']

aa_dfs = []
for df in dfs:
    aa_transporter = df[df.index.isin(elchain_genes)] #aa_list 
    aa_dfs.append(aa_transporter)

for x in aa_dfs:
    if '1_proCTC_rest' in x['label'].values:
        df= x.query('comparison == "g1_vs_g0"')
        print(df)

def clean_name(name):
    return re.sub(r'[\\/*?:\[\]]', '_', str(name))[:31]

with pd.ExcelWriter(os.path.join(path_results,'el_chain_proCTC_DEgs.xlsx'), engine='openpyxl') as writer:  #aa_transporters_in_DEgs.xlsx
    for df in aa_dfs:
        label_value= df['label'].iloc[0]
        sheet_name = clean_name(label_value)

        df_fil = df[['rank', 'evidence','effect_size','comparison','label']].copy()
        df_fil.rename(columns={'effect_size':'log2FC'}, inplace= True)

        df_fil.to_excel(writer, sheet_name=sheet_name)




#Search in Deg PPARGC1A and prepare excel contrast CTC vs lung
gene = ['PPARGC1A','PGC1','PGC1A','PGC-1alpha','PGC1alpha']
print(dfs)
for df in dfs:
    if 'label' in df.columns and (df['label'] == 'lung_ctc').any():
        lung_df = df[df['label'] == 'lung_ctc']
        filtered_df = lung_df[lung_df['comparison'] == 'g1_vs_g0']
        break

filtered_df

filtered_df = filtered_df[['evidence', 'effect_size']]
filtered_df = filtered_df.rename(columns={
    'evidence': 'FDR',
    'effect_size': 'logFC'
})

filtered_df.index.name = 'genes'
filtered_df = filtered_df.sort_values(by='logFC', ascending=False)
filtered_df.to_excel(os.path.join(path_results, 'DEGs_CTC_vs_Lung.xlsx'))

subset_df = filtered_df.loc[filtered_df.index.intersection(gene)]
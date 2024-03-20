"""
Manual annotation of cell states.
"""

import os
import pickle
from Cellula.plotting._colors import *
from Cellula.plotting._plotting_base import *
from Cellula.plotting._plotting import *
from BC_chemo_utils.tests import *
matplotlib.use('macOSX')


##


# Paths
path_main = '/Users/IEO5505/Desktop/BC_chemo_reproducibility'
path_data = os.path.join(path_main, 'data', 'MDA')
path_results = os.path.join(path_main, 'results', 'final_cell_states')


##


# Data
adata = sc.read(os.path.join(path_data, 'clustered.h5ad'))
embs = pd.read_csv(os.path.join(path_data, 'full_embs.csv'), index_col=0)
cell_state_map = (
        pd.read_csv(
        os.path.join(path_data, 'cell_states_definitive.csv'), 
        index_col=0
    )['Cell state'].to_dict()
)


##


# Here we go

# Check markers and paths
# c = 18
# comparison = f'{c}_vs_rest'
# 
# df['tmp'] = np.where(df.leiden == c, str(c), 'other')
# df['tmp'].value_counts()
# df['tmp'].value_counts(normalize=True)
# 
# d['df'].query('comparison == @comparison').head(10)
# d['gs'][comparison].compute_GSEA()
# d['gs'][comparison].GSEA['original'].sort_values('NES', ascending=False).head()
# 
# df.loc[:,conts+['tmp']].groupby('tmp').agg('median')
# 
# cov = 'condition'
# for x in df[cov].unique():
#     print(x)
#     compute_enrichment(df, 'leiden', cov, x).loc[c].to_dict()
#     print('\n')


##


d_annot = {

    0 : '''
    6904 cells. .12 entire dataset. 
    Stable to bootstrapping (>.75 support). At least 3 good marker genes.
    Proposed label: Cycling.
    Major cell cycle markers (CDK1, TOP2A, PCLAF). GSEA enrichment for 'DNA replication (GO:0006260)'
    and 'sister chromatid segregation pathways (GO:0000819)'. Coherent with cc signatures scores 
    (cycling 0.58 vs 0.05 rest, other cc signatures accordingly).
    No differences in nUMIs, mito_perc, apoptosis or n of detected genes with the rest 
    of the cells.
    Significantly enriched in NT_NT condition (.67 cells, FDR 0.0, OR 4.34)
    and PTs (.69 cells, FDR <10^-147, OR 1.98).
    ''',

    1 : '''
    6501 cells. .11 entire dataset. 
    Unstable to bootstrapping (<.60 support). No good markers. 
    Proposed label: IFN/Antigen presentation/TGFbeta.
    Shares IFI27, ISG15 and FN1 with other clusters. Other markers are PMEPA1 and SERPINF1
    (two interesting genes), SOX4, KRT7, LY6E. GSEA enrichment for type I interferon signaling
    pathway (GO:0060337). Similar (PAGA) to cluster 7, 5 and 14 (also, unstable and poorly
    characterized clusters). cc arrested (cycling 0.008393 vs 0.147665 rest, other cc 
    signatures accordingly).
    No differences in nUMIs, mito_perc, apoptosis or n of detected genes with the rest 
    of the cells.
    Significantly enriched in NT_NT and NT_AC conditions (.36 and .26 cells, 
    respectively, FDR <10^-16 and 10^-63, respectively, OR 1.29 and 1.69, respectively)
    and PTs (.93 cells, FDR~0, OR 13.72).
    ''',

    2: '''
    6211 cells. .10 entire dataset. 
    Least stable partition (~.45 support). No good markers. 
    Proposed label: Ribosome biogenesis
    Shares C12orf75 with other clusters. Other markers are AKAP12, PTTG1, FABP5,
    (interesting gene, FA binding protein) and CDKN3.
    GSEA enrichment for 'ribosome biogenesis (GO:0042254)' pathway (not so consistent with
    ribosomal gene signature upregulation, 2.96 vs 2.92 rest).
    Similar (PAGA) to cluster 10 and 14 (also, unstable and poorly characterized clusters). 
    Active cell cycle (cycling 0.24 vs 0.07 rest, other cc signatures accordingly).
    No major differences in nUMIs, mito_perc, apoptosis or n of detected genes with the rest 
    of the cells.
    Significantly enriched in NT_NT condition (.66 cells, FDR ~0, OR 3.97).
    Enriched both at PT and lungs (.62 cells, FDR <10^-30, OR 1.36 and 
    .38 cells, FDR <10^-07, OR 1.15, respectively).
    ''',

    3: '''
    5834 cells. .10 entire dataset. 
    Least stable partition (~.58 support). Not so good markers. 
    Proposed label: Endocytosis/MT-translation
    Shares PTTG1 with cluster 2. Other Marker genes UBE2C, UBE2S, AP2S1, CKS1B.
    GSEA enrichment for mitochondrial translation (GO:0032543).
    Similar (PAGA) to cluster 16 and 2. 
    Active cell cycle (cycling 0.36 vs 0.06 rest, other cc signatures accordingly).
    No major differences in nUMIs, mito_perc, apoptosis or n of detected genes with the rest 
    of the cells.
    Significantly enriched in NT_AC and AC_AC conditions (.45 and .32 cells, 
    respectively, FDR~0 and <10^-106, respectively, OR 4.42 and 1.97, respectively)
    and PTs (.93 cells, FDR~0, OR 13.72).
    Extremely enriched at lungs!! (.91 cells, FDR~0, OR 26.7 and 
    .38 cells, FDR <10^-07, OR 1.15, respectively).
    Note: if the question is, what cell state reaches metastasis after treatment, 
    here it might be an initial clue.
    ''',


    4: '''
    4593 cells. .08 entire dataset. 
    Stable partition (~.58 support). Relatively good markers.
    Proposed label: Ras/Survival
    Top markers: PAEP, CD59 (interesting gene). Other markers are TMEM158, PGK1 and SRPX
    (very interesting genes).
    All negative NES (NB: only GO Biological Process 2021).
    CC arrested (cycling 0.04 vs 0.11 rest, other cc signatures accordingly).
    No major differences in nUMIs, mito_perc, apoptosis or n of detected genes with the rest 
    of the cells.
    Significantly enriched in AC_AC condition (.54 cells, FDR~0, OR 5.46)
    and PTs (.93 cells, FDR~0, OR 13.72).
    Extremely enriched in PTs (.86 cells, FDR~0, OR 5.26).
    Note: Enriched in PTs but not in AC_NT condition... Why?
    ''',

    5: '''
    4180 cells. .07 entire dataset. 
    Stable partition (~.63 support). Not so good markers.
    Proposed label: IFN/undefined
    Shares MALAT1, NEAT1, SOX4 NNMT and SERPINF1 with other clusters.
    GSEA enrichment for type I interferon signaling pathway (GO:0060337).
    CC arrested (cycling -0.01 vs 0.14 rest, other cc signatures accordingly).
    No major differences in mito_perc, apoptosis or n of detected genes with the rest 
    of the cells. 13k UMIs vs 21k rest (median).
    Significantly enriched in AC_AC condition (.36 cells, FDR~0, OR 3.71)
    and PTs (.87 cells, FDR~0, OR 5.69).
    Note: As before. Cluster 5 Enriched in AC_NT PTs but not in AC_AC ones... Why?
    ''',

    6: '''
    3064 cells. .05 entire dataset. 
    Stable partition (>~.8 support). Very good markers: CST1,2,4.
    Proposed label: Cystatine/Ca2+-metabolism
    GSEA enrichment for 'negative regulation of calcium ion transmembrane transporter 
    activity (GO:1901020)'.
    Active cell cycle entry, but no G2/M checkpoint passing (cycling 0.11 vs 0.10 rest,
    cycle_diff 0.04 vs 0.06, s_seurat 0.005 vs -0.01 rest, g2m_seurat 0.036 vs 0.044 rest).
    No major differences in mito_perc, apoptosis or n of detected genes with the rest 
    of the cells.
    Significantly enriched in both AC_AC condition (.4 cells, FDR~0, OR 2.61) and NT_AC
    (.25 cells, FDR<10^-17, OR 1.45) and lungs (.69 cells, FDR~0, OR 4.21). 
    Very interesting cluster.
    ''',

    7: '''
    2589 cells. .04 entire dataset. 
    Stable partition (>~.75 support). Very good markers.
    Proposed label: Immune-like
    Top markers are CEBPD, NNMT, CXCL8, CXCL2, JUNB, partially shared with cluster 5. 
    SOX4 shared with other clusters.
    GSEA enrichment for: 
    - 'cellular response to cytokine stimulus (GO:0071345)',
    - 'cellular response to type I interferon (GO:0071357)',
    - 'neutrophil chemotaxis (GO:0030593)',
    - 'neutrophil migration (GO:1990266)'.
    CC arrest (cycling -0.04 vs 0.12 rest).
    Difference in nUMIs 14k vs 21 k (median) others, but not in other covariates.
    Significantly enriched in NT_NT condition (.56 cells, FDR~0, OR 2.1) and PT 
    (.68 cells, FDR<10^-26, OR 1.56).
    Note: Same question as some months ago. Is a neutrophil program aberrantly activated
    by tumor cells, or a tumor response to neutrophils?
    ''',

    8: '''
    2535 cells. .04 entire dataset. 
    Stable partition (>~.75 support). Very good markers
    Proposed label: IFN/TNFalpha
    CHE DIO MI STRAFULMINI, i top3 markers sono:
    NFKBIA, ANGPTL4 e IFI27. ANGPTL4 con 1.740663 log2FC, .94 perc_group
    e .54 perc rest. Other relevant markers: TNFAIP3. NFKBIA and CXCL2 are shared 
    with cluster 7.
    GSEA enrichment for interferon-gamma-mediated signaling pathway (GO:0060333) and other IFN
    pathways.
    CC arrest (cycling -0.01 vs 0.12 rest). No major differences in mito_perc, 
    apoptosis or n of detected genes with the rest of the cells.
    Significantly enriched in AC_NT condition (.95 cells, FDR~0, OR 157.22) and PT 
    (.99 cells, FDR~0, OR 76.39).
    ''',

    9: '''
    2446 cells. .04 entire dataset. 
    Not so stable partition (<.7 support). Good markers.
    Proposed label: AP1
    Top markers FOS, JUN, EGR1. Other GADD45A,B IER1,5 (interesting, both), partially shared
    with 5,7 and most of all 18 clusters.
    GSEA enrichment are all negative.
    Active cc and cell division (cycling 0.33 vs 0.05 rest, all other signatures accordingly). 
    No major differences in mito_perc, apoptosis or n of detected genes with the rest of the cells.
    Significantly enriched in NT_NT condition (.57 cells, FDR~0, OR 2.11), 
    NT_AC condition (.25 cells, FDR<10^-10, OR 1.37) and lungs (.87 cells, FDR~0, OR 13.22).
    ''',

    10: '''
    2329 cells. .04 entire dataset. 
    Not so stable partition (<.6 support). Poor markers.
    Proposed label: Ribosomecbiogenesis
    Only one marker gene, SNHG25. Others are shared.
    GSEA enrichment for: 'ribosome biogenesis (GO:0042254)', 'rRNA metabolic process (GO:0016072)'.
    Active cc and cell division (cycling 0.22 vs 0.05 rest, all other signatures accordingly). 
    No major differences in mito_perc, apoptosis or n of detected genes with the rest of the cells.
    Significantly enriched in NT_AC condition (.63 cells, FDR~0, OR 7.72) 
    and lungs (.92 cells, FDR~0, OR 22.3).
    ''',

    11: '''
    2180 cells. .04 entire dataset. 
    Stable partition (>.8 support). Good markers.
    Proposed label: IFN
    Very related to cluster 12. IGFBP6 (very interesting gene) is the only specific genes 
    with respect to cluster 12. Others (ISG15,20, OASL, LY6E, PLAAT4, IFIT3) 
    are shared at some level.
    GSEA enrichment for 'cellular response to type I interferon (GO:0071357)'.
    CC arrested (cycling 0.02 vs 0.11 rest, all other signatures accordingly). 
    Also the s_seurat signature seems pretty downregulated.
    No major differences in mito_perc, apoptosis or n of detected genes with the rest of the cells.
    Significantly enriched in AC_AC (.24 cells, FDR~0.01, OR 1.18) and AC_NT 
    (.24 cells, FDR<10^-29, OR 1.83) conditions, and lungs (.44 cells, FDR<10^-11, OR 1.34).
    Note: it might some residual IFN state that remains activated after metastatization?
    Is it truly a state that emerges only after PT treatment?? Quiescence must be verified.
    ''',

    12: '''
    2116 cells. .04 entire dataset. 
    Stable partition (>.8 support). Very good markers.
    Proposed label: IFN/growth factor response.
    All IFIT genes (2,3,1), OASL and ISG15 (shared with cluster 11).
    Other relevant genes: KLF4 (Yamanaka pluripotency factor), ZC3HAV1 and ZFP36L2. 
    GSEA enrichment for 'defense response to symbiont (GO:0140546)' and 
    'cellular response to type I interferon (GO:0071357)'.
    CC arrested (cycling 0.02 vs 0.11 rest, all other signatures accordingly). 
    Also the s_seurat signature seems pretty downregulated.
    No major differences in mito_perc, apoptosis or n of detected genes with the rest of the cells.
    Significantly enriched in AC_AC (.24 cells, FDR~0.02, OR 1.13) and AC_NT 
    (.24 cells, FDR<10^-29, OR 1.83) conditions, and PTs (.74 cells, FDR<10^-54, OR 2.11).
    Note: How much cluster 12 and 12 are coupled? Are they the states in which the same clones
    evolve in this two branch across the PT-lung trajectory? How much are they similar??
    What are their differences?
    ''',

    13: '''
    1757 cells. .03 entire dataset. 
    Stable partition (>.8 support). Very good markers.
    Proposed label: Glicolysis/Stress-response/Ras-mediated-survival
    A lot of interesting genes: SLC2A1, NDRG1, ADM, ENO2 and TMEM158.
    GSEA enrichment for 'carbohydrate catabolic process (GO:0016052)' and 
    'glucose catabolic process to pyruvate (GO:0061718)'.
    CC arrested (cycling -0.01 vs 0.12 rest, all other signatures accordingly). 
    Also the s_seurat signature seems pretty downregulated.
    No major differences in mito_perc, apoptosis or n of detected genes with the rest 
    of the cells. Low nUMIs (15k vs 21k rest) and detected genes (3.9k vs 4.9k rest).
    Significantly enriched in AC_AC (.29 cells, FDR<10^-14, OR 1.51) and AC_NT 
    (.39 cells, FDR<10^-129, OR 3.67) conditions, and PTs (.98 cells, FDR~0, OR 44.41).
    Note: It might be the glycolytic state we have found earlier.
    ''',

    14: '''
    1292 cells. .02 entire dataset. 
    Very unstable partition (<.6 support). Not so good markers.
    Proposed label: TGFBeta/Epitelial/Wound healing-like
    NEAT1, MALAT1, FN1, CEBPD, SOX4 markers shared with other clusters.
    ITGB4 within top10 markers. Other interesting gene: TGFB2.
    GSEA enrichment for 'cell junction assembly (GO:0034329)' and 
    'wound healing (GO:0042060)'.
    CC arrested (cycling -0.01 vs 0.11 rest, all other signatures accordingly). 
    Also the s_seurat signature seems pretty downregulated.
    No major differences in mito_perc, apoptosis or n of detected genes with the rest 
    of the cells. Low nUMIs (16k vs 21k rest).
    Significantly enriched in NT_NT (.9 cells, FDR~04, OR 14.62)
    condition, and lungs (.61 cells, FDR~0, OR 2.65).
    Note: Might be cells that already underwent MET after the metastatic cascade??
    ''',

    15: '''
    865 cells. .01 entire dataset. 
    Very stable partition (>.8 support). Good markers.
    Proposed label: Ribosome biogenesis/FA metabolism/GAS5-S100A4
    GAS5 is the top marker! Other interesting gene: S100A4, RHOBTB3 and BST2. 
    Not in top10 markers, but highly specific and potentially relevant:
    TFF3 and DMKN genes.
    GSEA enrichment for 'ribosome biogenesis (GO:0042254)' and 
    'acylglycerol metabolic process (GO:0006639)'.
    CC arrested (cycling -0.01 vs 0.11 rest, all other signatures accordingly). 
    Also the s_seurat signature seems pretty downregulated.
    No major differences in mito_perc, apoptosis or n of detected genes with the rest 
    of the cells. High nUMIs (30k vs 20k rest).
    Significantly enriched in AC_AC (.74 cells, FDR~0, OR 10.52)
    condition, and PTs (.64 cells, FDR~0, OR 64.19).
    ''',

    16: '''
    816 cells. .01 entire dataset. 
    Very stable partition (>.8 support). Relatively good markers.
    Proposed label: OXPHOS/C12orf75
    Interesting genes C12orf75, PTTG1.
    GSEA enrichment for 'aerobic electron transport chain (GO:0019646)'.
    Actively cycling (cycling 0.15 vs 0.10 rest, all other signatures accordingly). 
    Very low nUMIs (13k vs 20k rest) and genes (2.7k vs 4.9k).
    Significantly enriched in NT_AC (.36 cells, FDR<10^-26, OR 2.26)
    condition, and lungs (.56 cells, FDR<10^-25, OR 2.1).
    ''',

    17: '''
    354 cells. .006 entire dataset. 
    Very stable partition (>.8 support). Best markers.
    Proposed label: Chemotaxis/Migration
    Interesting genes ZNF90, a lot of H3F3 histon variants.
    GSEA enrichment for:
    -'positive chemotaxis (GO:0050918)'
    -'TRIF-dependent toll-like receptor signaling pathway (GO:0035666)'
    -'MyD88-dependent toll-like receptor signaling pathway (GO:0002755)'
    CC arrested (cycling 0.02 vs 0.10 rest, all other signatures accordingly). 
    Very low nUMIs (3k vs 20k rest) and genes (<1k vs 4.8k rest).
    Significantly enriched in AC_NT (.92 cells, FDR~0, OR 66.08)
    condition, and lungs (.97 cells, FDR~0, OR 55.68).
    Note: This cluster is the shittiest for QC metrics, and is enriched for the shitty branch.
    ''',

    18: '''
    227 cells. .003 entire dataset. 
    Most stable partition (>.8 support). Very good markers.
    Proposed label: TGFBeta/ECM
    Interesting genes SPARC, COL3A1, ID3. Shares FOS, JUN, EGR1 with cluster 9. 
    GSEA enrichment for:
    -'cellular response to transforming growth factor beta stimulus (GO:0071560)'
    -'positive regulation of pri-miRNA transcription by RNA polymerase II (GO:1902895)'
    Actively cycling (cycling 0.14 vs 0.10 rest, all other signatures accordingly). 
    No major QC metrics differences.
    Significantly enriched in NT_NT (.95cells, FDR~0.003, OR 1.51)
    condition, and lungs (.83 cells, FDR~0, OR 7.95).
    '''

}

# Save
# (
#     pd.Series(d_annot, name='annotation')
#     .reset_index()
#     .rename(columns={'index':'cluster'})
#     .to_csv(os.path.join(path_results, 'annotation.csv'), index=False)
# )


##


# # Write df.obs cols
# manual_annot = {
# 
#     0 : 'Cycling',
#     1 : 'IFN/Antigen presentation/TGFbeta',
#     2: 'Ribosome biogenesis',
#     3: 'Endocytosis/MT-translation',
#     4: 'Ras/Survival',
#     5: 'IFN/undefined',
#     6: 'Cystatine/Ca2+-metabolism',
#     7: 'Immune-like',
#     8: 'IFN/TNFalpha',
#     9: 'AP1',
#     10: 'Ribosome biogenesis',
#     11: 'IFN',
#     12: 'IFN/growth factor response',
#     13: 'Glicolysis/Stress-response/Ras-mediated-survival',
#     14: 'TGFBeta/Epitelial/Wound healing-like',
#     15: 'Ribosome biogenesis/FA metabolism/GAS5-S100A4',
#     16: 'OXPHOS/C12orf75',
#     17: 'Chemotaxis/Migration',
#     18: 'TGFBeta/ECM'
# 
# }

# manual_annot_small = {
# 
#     0 : 'Cycling',
#     1 : 'IFN',
#     2: 'Ribosome biogenesis',
#     3: 'Endocytosis',
#     4: 'Ras/Survival',
#     5: 'IFN',
#     6: 'Cystatine/Ca2+-metabolism',
#     7: 'Immune-like',
#     8: 'IFN',
#     9: 'AP1',
#     10: 'Ribosome biogenesis',
#     11: 'IFN',
#     12: 'IFN',
#     13: 'Glicolysis',
#     14: 'TGFBeta',
#     15: 'Ribosome biogenesis',
#     16: 'OXPHOS',
#     17: 'Chemotaxis/Migration',
#     18: 'TGFBeta'
# 
# }

# Re-annotate
# df['cell_state_highres'] = df['leiden'].map(lambda x: manual_annot[x])
# df['cell_state_lowres'] = df['leiden'].map(lambda x: manual_annot_small[x])


# New, final annot
adata.obs['leiden'] = adata.obs['leiden'].astype('int')
adata.obs['final_cell_state'] = adata.obs['leiden'].map(cell_state_map)

# Save clustered and full_embs
adata.write(os.path.join(path_data, 'clustered.h5ad'))
embs['final_cell_state'] = adata.obs['final_cell_state']
embs.to_csv(os.path.join(path_data, 'full_embs.csv'))


##
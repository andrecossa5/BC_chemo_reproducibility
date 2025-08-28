import pandas as pd
import gseapy as gp
import gc
import os
import pandas as pd 
from gseapy import prerank
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
from plotting_utils.plotting_base import *
matplotlib.use('macOSX')

class GSEAAnalysis:
    def __init__(self, pca_loadings_df, organism='human'):
        """
        Initialize the GSEA Analysis class with PCA loadings.
        """
        self.pca_loadings_df = pca_loadings_df  
        self.organism = organism
        self.GSEA = {}

    def compute_GSEA(self, covariate='PC1_loading', collection='GO_Biological_Process_2023'):
        """
        Perform GSEA using PCA loadings (e.g., PC1 or PC2) as covariate.
        """
        if covariate not in self.pca_loadings_df.columns:
            raise ValueError(f"{covariate} not found in the PCA loadings dataframe")

        self.pca_loadings_df['gene'] = self.pca_loadings_df['gene'].astype(str)

        self.pca_loadings_df = self.pca_loadings_df.dropna(subset=[covariate])

        ranked_gene_list = self.pca_loadings_df[['gene', covariate]].sort_values(by=covariate, ascending=False)
        ranked_gene_list.set_index('gene', inplace=True)

        results = prerank(
            rnk=ranked_gene_list,
            gene_sets=[collection],
            threads=4,
            min_size=15,
            max_size=500,
            permutation_num=200, 
            outdir=None, 
            seed=1234,
            verbose=True,
        )

        df = results.res2d.loc[:, 
            [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
        ].rename(columns={'FDR q-val' : 'Adjusted P-value'})
        
        idx = df.sort_values(by='Adjusted P-value').index
        filtered_df = df.loc[idx, :]

        filtered_df['Term'] = filtered_df['Term'].map(lambda x: str(x).split('__')[1] if isinstance(x, str) else str(x))
        filtered_df = filtered_df.set_index('Term')

        self.GSEA['original'] = filtered_df
        gc.collect()
        
        return filtered_df

path="/Users/ieo7295/Desktop/BC_sh/data"
pca_loadings_df = os.path.join(path,"pca_loadings.csv")
pca_loadings_df=pd.read
gsea_analysis = GSEAAnalysis(pca_loadings_df)

#Run the GSEA
gsea_results_pc1 = gsea_analysis.compute_GSEA(covariate='PC1_loading', collection='GO_Biological_Process_2023')
gsea_results_pc2 = gsea_analysis.compute_GSEA(covariate='PC2_loading', collection='GO_Biological_Process_2023')

gsea_results_pc1.to_excel('GSEA_PC1.xlsx')
gsea_results_pc2.to_excel('GSEA_PC2.xlsx')


"""
Plot stem-plots of top 2 PCs GSEA-enriched pathways
"""
i=1 
pc_column = f'PC{i}_loading'  
loading_data = pca_loadings_df[['gene', pc_column]].sort_values(by=pc_column, ascending=False)
organism='human'
collection='GO_Biological_Process_2023'
gsea = GSEAAnalysis(loading_data, organism=organism)
gsea_results = gsea.compute_GSEA(covariate=pc_column, collection=collection)

   
fig, axs = plt.subplots(1, 2, figsize=(15, 8))
stem_plot(
        gsea_results.iloc[:, [0, 1, 3]].sort_values('NES', ascending=False).head(50),
        'NES', 
        ax=axs[0]
    )
format_ax(axs[0], title='GSEA', xlabel='NES')
  
stem_plot(
        loading_data[f'PC{i}_loading'].to_frame('es').sort_values('es', ascending=False).head(50),
        'es', 
        ax=axs[1]
    )
format_ax(axs[1], title=f'PC{i} loadings', xlabel='Loadings')
fig.suptitle(f'PC{i}, scaled layer, original representation')
fig.tight_layout()
#fig.savefig((".png"),dpi=300)
plt.show()

gsea_results=df_deg

### GSEA for regulatory network genes ###

def gene_gsea(df_deg ,collection='GO_Biological_Process_2023'):

     df_deg['gene'] = df_deg['gene'].astype(str)

     ranked_gene_list = df_deg.sort_values(by='logFC', ascending=False)
     ranked_gene_list = ranked_gene_list.set_index('gene')['logFC']

     results = prerank(
            rnk=ranked_gene_list,
            gene_sets=[collection],
            threads=4,
            min_size=15,
            max_size=500,
            permutation_num=200, 
            outdir=None, 
            seed=1234,
            verbose=True,
        )

     df = results.res2d.loc[:, 
            [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
        ].rename(columns={'FDR q-val' : 'Adjusted P-value'})
    
     idx = df.sort_values(by='Adjusted P-value').index
     filtered_df = df.loc[idx, :]

     filtered_df['Term'] = filtered_df['Term'].map(lambda x: str(x).split('__')[1] if isinstance(x, str) else str(x))
     filtered_df = filtered_df.set_index('Term')
     gc.collect()
        
     return filtered_df


file_path="/Users/ieo7295/Desktop/BC_chemo_reproducibility"
alt_file_path="/Users/ieo7295/Desktop/BC_sh/results/res_final"
file=os.path.join(file_path,'results', 'MDA','Gsea_promet_AC_vs_promet_NT.csv')
alt_file=os.path.join(alt_file_path,"Degs_scr_vs_paep_rm_out_new.xlsx")
df_deg=pd.read_excel(alt_file)
df_deg= df_deg.rename(columns={'Unnamed: 0':'gene'})
df_deg=df_deg.drop(columns={'Unnamed: 0'})
df_deg_pvalue=pd.read_excel(alt_file)
df_deg_pvalue=df_deg_pvalue.rename(columns={'Unnamed: 0':'gene'})

network_gene_2=pd.read_csv("PAEP1vsSCR_gsea.csv")
network_gene_2= network_gene_2.rename(columns={'Unnamed: 0':'gene'})

network_gene_3=pd.read_csv("PAEP2vsSCR_gsea.csv")
network_gene_3= network_gene_3.rename(columns={'Unnamed: 0':'gene'})

gsea = gene_gsea(df_deg,collection='GO_Biological_Process_2023')
gsea.to_excel(os.path.join(alt_file_path,"SCRvsshPAEP_gsea_new_15_500.xlsx"))

gsea.to_excel("PAEP1vsSCR_gsea.xlsx", index=True)
gsea.to_excel("PAEP2vsSCR_gsea.xlsx", index=True)
gsea.head(50)

gsea = gsea.query('`Adjusted P-value` <= 0.1')
fig, ax = plt.subplots(figsize=(15, 8))
stem_plot(
    gsea.iloc[:, [0, 1, 3]].head(50),
    'NES', 
    ax=ax 
)
format_ax(ax, title='GSEA NES', xlabel='NES')
ax.set_xlim(-2.5,2.5)
ax.set_xticks(np.arange(-2.5,2.6,0.5))
plt.tight_layout()
plt.show()
plt.savefig(os.path.join(alt_file_path, "GSEA_NES_scr_paep_new.png"),dpi=300)


### gene enrichment plot ###
df_pc1=pd.read_excel('GSEA_PC1.xlsx')
df_pc2=pd.read_excel('GSEA_PC2.xlsx')


def plot_gsea_results(df_deg, title='Gene Set Enrichment Analysis'):
    #df = pd.read_excel(df_deg, index_col=0)
    df=df_deg
    df = df.reindex(df['NES'].head(50))
    
    plt.figure(figsize=(10, 6))
    sns.barplot(y=df.index, x=df['NES'], palette='viridis', edgecolor='black')
    
    plt.axvline(x=0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.ylabel('Gene Sets')
    plt.title(title)
    
    plt.tight_layout()
    plt.show()

# Example usage
plot_gsea_results(df_deg, title='GSEA Results for PC1 Loading')
df_deg= df_deg.set_index('Term')
def plot_GSEA_dot(gsea, title="Enrichment Dot Plot", top_n=50):
    gsea['Adjusted P-value'] = pd.to_numeric(gsea['Adjusted P-value'], errors='coerce')
    df_plot = gsea.nsmallest(top_n, 'Adjusted P-value')

    plt.figure(figsize=(14, 9))
    scatter = plt.scatter(
        x=df_plot['Adjusted P-value'],
        y=df_plot.index,  
        c=df_plot['Adjusted P-value'],
        cmap='coolwarm',
        edgecolors='black'
    )
    
    plt.xlabel('Adjusted P-value')
    plt.ylabel('Pathway')
    plt.colorbar(scatter, label="Adjusted P-value")
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(file_path, "GSEA_deg_paepvsscr_dotplot_new.png"),dpi=300)
    plt.show()

plot_GSEA_dot(df_deg)
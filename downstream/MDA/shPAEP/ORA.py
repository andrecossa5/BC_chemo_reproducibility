import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gc
import os
import networkx as nx
from gseapy import enrichr, prerank

class ORAAnalysis:
    def __init__(self, subset_deg, organism='human'):
        self.subset_deg= subset_deg
        self.organism = organism
        self.ORA_results = {}


    def compute_ORA(self, key='ORA_results', gene_column='Gene',
        collection='GO_Biological_Process_2023', n_out=100):
        """
        Perform ORA (Over-Representation Analysis)
        """
        if gene_column not in self.subset_deg.columns:
            raise ValueError(f"Column '{gene_column}' not found in dataframe")
        
        gene_list = self.subset_deg[gene_column].dropna().tolist()

        results = enrichr(
            gene_list=gene_list,
            gene_sets=[collection],
            organism=self.organism, 
            outdir=None, 

        ).results
        results.to_excel(os.path.join(file_path, "ORA_results_deg_regulon_genes.xlsx"))
        df = results.loc[:, 
            [ 'Term','Overlap','Adjusted P-value', 'Genes' ]
        ]

  
        filtered_df = df.nsmallest(n_out, 'Adjusted P-value')
        filtered_df = filtered_df.set_index('Term')

 
        self.ORA_results[key] = filtered_df

 
        gc.collect()

        return filtered_df 

file_path="/Users/ieo7295/Desktop/BC_sh/results/res_final"
file=os.path.join(file_path,"Degs_paep_vs_scr_rm_out_thrb_regulon.xlsx")
file_genes=os.path.join(file_path,"allgenes_paepvsscr.xlsx")
df_deg=pd.read_excel(file)
df_genes=pd.read_excel(file_genes)
df_deg.rename(columns={'Unnamed: 0':'Gene'}, inplace=True)
df_genes.rename(columns={'Unnamed: 0':'Gene'}, inplace=True)
#subset_deg= pd.concat([df_deg.head(50), df_deg.tail(50)])
subset_deg=df_deg
ora_analysis= ORAAnalysis(subset_deg)


ora_results = ora_analysis.compute_ORA(gene_column='Gene')
ora_results.to_excel(os.path.join(file_path, "ORA_results_thrb_regulon_shPAEPvsSCR_hallmark.xlsx"))


def plot_ORA_bar(ora_results, title="Top Enriched Pathways", top_n=50):
    """
    Plots a bar chart of the top enriched pathways based on Adjusted P-value.

    """
    # Select top pathways
    df_plot = ora_results.nsmallest(top_n, 'Adjusted P-value')

    plt.figure(figsize=(13, 9))
    sns.barplot(
        data=df_plot,
        y=df_plot.index,
        x='Adjusted P-value',
        palette='viridis'
    )
    
    plt.xlabel('Adjusted P-value', fontsize=16)
    plt.ylabel('Pathway',)
    plt.xlim(left=0.04)
    plt.title(title) 
    plt.tight_layout()
    plt.savefig(os.path.join(file_path, "ORA_deg_paep_vs_scr_rm_out_thrbregulon_plotbar.png"),dpi=300)



plot_ORA_bar(ora_results)



def plot_ORA_dot(ora_results, title="Enrichment Dot Plot", top_n=50):
    """
    Plots a dot plot of enriched pathways showing adjusted P-value and gene overlap.
    
    """

    df_plot = ora_results.nsmallest(top_n, 'Adjusted P-value')
    
    plt.figure(figsize=(15, 10))
    scatter = plt.scatter(
        x=df_plot['Adjusted P-value'],
        y=df_plot.index,
        s=df_plot['Overlap'].apply(lambda x: int(x.split('/')[0])) * 50,  # Bubble size = overlapping genes
        c=df_plot['Adjusted P-value'],
        cmap='coolwarm',
        edgecolors='black'
    )
    
    plt.xlabel('Adjusted P-value',fontsize=16)
    plt.ylabel('Pathway',fontsize=16)
    plt.colorbar(scatter, label="Adjusted P-value")
    plt.title(title,fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(file_path, "ORA_deg_paep_vs_scr_rm_out_thrbregulon_dotplot_hallmark.png"),dpi=300)



plot_ORA_dot(ora_results)





def plot_enrichment_map(ora_results, top_n=50):
    """
    Creates a network graph of enriched pathways.
    """
    df_plot = ora_results.nsmallest(top_n, 'Adjusted P-value')

    G = nx.Graph()
    
    for _, row in df_plot.iterrows():
        pathway = row.name  
        genes = row['Genes'].split(";")  
        
        G.add_node(pathway, type='pathway')
        
        for gene in genes:
            G.add_node(gene, type='gene')
            G.add_edge(pathway, gene)
    
    # Draw the network
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    node_colors = ['red' if G.nodes[n]['type'] == 'pathway' else 'blue' for n in G.nodes()]
    
    nx.draw(
        G, pos, with_labels=True, node_color=node_colors, edge_color='gray',
        node_size=500, font_size=8, alpha=0.8
    )
    plt.title("Enrichment Map")
    plt.show()


plot_enrichment_map(ora_results)

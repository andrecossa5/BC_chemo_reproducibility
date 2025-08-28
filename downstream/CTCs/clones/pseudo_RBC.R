# DE with edgeR
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)
library(edgeR)
library(readxl)
library(fgsea)
library(AUCell)
library(GSVA)

##

#Utils: aggregate data into pseudobulk per group
aggregate_data <- function(sce, assay='raw', group_by=NULL, cont_covariates=c('nUMIs', 'mito_perc'), 
                           cat_covariates=c('seq_run')) {
  
  meta <- colData(sce) %>% as.data.frame()
  cats <- levels(meta[[group_by]])
  
  M <- assays(sce)[[assay]]
  L <- lapply(cats, function(x) { rowSums( M[, meta[[group_by]] == x] %>% as.matrix() ) })
  L <- setNames(L, cats)
  df <- as.data.frame(L) %>% t()
  colnames(df) <- row.names(sce)
  
  other_cont <- meta %>% 
    select(all_of(c(group_by, cont_covariates))) %>% 
    group_by(.data[[group_by]]) %>% 
    summarise_all(median, na.rm=TRUE)
  other_cat <- meta %>% 
    select(all_of(c(group_by, cat_covariates))) %>% 
    group_by(.data[[group_by]]) %>%
    unique() 
  
  df_cov <- merge(other_cont, other_cat, by=group_by)
  df <- merge(df_cov, df %>% as.data.frame() %>% rownames_to_column(var=group_by), by=group_by)
  row.names(df) <- df[[group_by]]
  df <- df %>% select(-c(group_by))
  
  return(df)
}

#Fit pseudobulk with edgeR fixed effect model
fit_pseudobulk_fixed <- function(df_pseudo, sce, test_column=NULL, design=NULL) {
  y <- DGEList(
    counts = t(df_pseudo[, colnames(df_pseudo) %in% row.names(sce)]), 
    samples = df_pseudo[, !colnames(df_pseudo) %in% row.names(sce)]
  )
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  y <- estimateDisp(y, design=design)
  fit <- glmQLFit(y, design)
  return(fit)
}

# Test contrasts with edgeR
test_pseudobulk_fixed <- function(fit, contrast=NULL) {
  qlf <- glmQLFTest(fit, contrast=contrast)
  tt <- topTags(qlf, n=Inf)$table
  tt <- tt %>% arrange(desc(logFC))
  return(tt)
}

### Here we go ###

# List of GBC barcodes of interest
gbc_list <- c("TATGGCGGTTGTGCTGTC", "CCATAAGGGGAGGATGTA", "TGACTGTGGAGTCCTGAA", 
              "TTTGTACGTGGGACCATA", "CGCGACTAGAGACGGCTC", "CTGCGGTTTCGTTAACGC", 
              "ATCTGGGACCTCAACAAG", "CGCGTCACACTGTCGGGC", "CTCGAGTCTGGATGCGTA", 
              "GGTGCAGGCGGTGTGGGG", "TTCCGAGCAGCTGCTTGG")
##

# Paths
path_main <- '/Users/ieo7295/Desktop/BC_chemo_reproducibility/'
path_data <- paste0(path_main, '/data/CTCs/grn/resources')
path_results <- paste0(path_main, 'results/CTCs/pseudobulk')

# Load Anndata as sce
sce <- readH5AD(paste0(path_data, "/clustered.h5ad"))


##                
                
                
colData(sce)$group <- paste0(colData(sce)$GBC, "_", colData(sce)$origin)

# Filter sce for GBC and origin in (CTC, PT, lung)
sce_filtered <- sce[, 
                    colData(sce)$GBC %in% gbc_list & 
                      colData(sce)$origin %in% c("CTC", "PT", "lung")
]

colData(sce_filtered)$group <- factor(colData(sce_filtered)$group)


min_cells <- 10
group_counts <- table(colData(sce_filtered)$group)
valid_groups <- names(group_counts[group_counts >= min_cells])
sce_filtered <- sce_filtered[, colData(sce_filtered)$group %in% valid_groups]
colData(sce_filtered)$group <- factor(colData(sce_filtered)$group)

# Aggregate to pseudobulk dataframe
df_pseudo <- aggregate_data(
  sce_filtered, 
  assay='raw', 
  group_by='group',
  cont_covariates=c('nUMIs', 'mito_perc','G1.S','G2.M'), 
  cat_covariates=c('seq_run')
)

df_pseudo <- df_pseudo %>%
  rownames_to_column("group") %>%
  separate(group, into=c("GBC", "origin"), sep="_", remove=FALSE)

df_pseudo$test_column <- factor(df_pseudo$origin, levels=c("CTC", "PT", "lung"))

# Design matrix
design <- model.matrix(~ 0 + test_column + nUMIs + mito_perc + G2.M, data=df_pseudo)

# Define contrasts for pro_PT_vs_CTC and pro_lung_vs_CTC
contrast_CTC_vs_PT <- makeContrasts(test_columnCTC - test_columnPT, levels=design)
contrast_CTC_vs_lung <- makeContrasts(test_columnCTC - test_columnlung, levels=design)
#contrast_CTC_vs_PT_lung <- makeContrasts(test_columnCTC - (test_columnPT + test_columnlung)/2, levels=design)
contrast_PT_vs_CTC <- makeContrasts(test_columnPT - test_columnCTC, levels=design)
contrast_lung_vs_CTC <- makeContrasts(test_columnlung - test_columnCTC, levels=design)

# Fit model
fit <- fit_pseudobulk_fixed(df_pseudo, sce_filtered, test_column=df_pseudo$test_column, design=design)

# Test contrasts
res_CTC_vs_PT <- test_pseudobulk_fixed(fit, contrast_CTC_vs_PT)
res_CTC_vs_lung <- test_pseudobulk_fixed(fit, contrast_CTC_vs_lung)
#res_CTC_vs_PT_lung <-test_pseudobulk_fixed(fit, contrast_CTC_vs_PT_lung)
res_PT_vs_CTC <- test_pseudobulk_fixed(fit, contrast_PT_vs_CTC)
res_lung_vs_CTC <- test_pseudobulk_fixed(fit, contrast_lung_vs_CTC)

# Save results (adjust paths as needed)
write.csv(res_CTC_vs_PT, paste0(path_results, '/results_pro_CTC_vs_PT.csv'), row.names=TRUE)
write.csv(res_CTC_vs_lung, paste0(path_results, '/results_pro_CTC_vs_lung.csv'), row.names=TRUE)
write.csv(res_CTC_vs_PT_lung, paste0(path_results, '/results_pro_CTC_vs_PT_lung.csv'), row.names=TRUE)
write.csv(res_PT_vs_CTC, paste0(path_results, '/results_pro_PT_vs_CTC.csv'), row.names=TRUE)
write.csv(res_lung_vs_CTC, paste0(path_results, '/results_pro_lung_vs_CTC.csv'), row.names=TRUE)
         
#Load RBC dataset 
bulk_data <- read_excel(paste0(path_main, '/data/CTCs/cell_states/GSE273783_Counts.xlsx'))
bulk_tpm <- read_excel(paste0(path_main, '/data/CTCs/cell_states/GSE273783_TPM.xlsx'))

bulk_expr <- bulk_tpm %>%
  filter(!is.na(Gene_Symbol)) %>%               
  distinct(Gene_Symbol, .keep_all = TRUE) %>%   
  column_to_rownames("Gene_Symbol") %>%        
  select(ends_with("_TPM")) %>%                
  select(-matches("CN|LIPO")) 

bulk_expr_raw <- bulk_data %>%
  filter(!is.na(Gene_Symbol)) %>%               
  distinct(Gene_Symbol, .keep_all = TRUE) %>%   
  column_to_rownames("Gene_Symbol") %>%        
  select(ends_with("_Count")) %>%                
  select(-matches("CN|LIPO")) 

group <- ifelse(grepl("CT", colnames(bulk_expr)), "healthy", "metastatic")

up_CTC <- res_lung_vs_CTC %>%
  filter(logFC > 1, FDR < 0.05) %>%
  rownames()

down_CTC <- res_lung_vs_CTC %>%
  filter(logFC < -1 , FDR < 0.05 ) %>%
  rownames()

gene_sets <- list(CTC_up = up_CTC, CTC_down = down_CTC)

#check overlap of genes
de_genes <- rownames(res_PT_vs_CTC)
bulk_genes <- bulk_data$Gene_Symbol
common_genes <- intersect(up_CTC, bulk_genes)
filtered_gene_sets <- lapply(gene_sets, function(genes) {
  intersect(genes, rownames(bulk_expr))
})


#DE bulk healthy vs metastatic
sample_names <- colnames(bulk_expr_raw)
group <- ifelse(grepl("MV", sample_names), "Metastatic", "Healthy")
group <- factor(group)
design <- model.matrix(~ 0 + group)

y <- DGEList(counts = bulk_expr_raw, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

contrast_metastatic_healthy <- makeContrasts(groupMetastatic - groupHealthy, levels=design)
res_met_vs_hlt <- test_pseudobulk_fixed(fit, contrast_metastatic_healthy)

contrast_healthy_metastatic <-  makeContrasts(groupHealthy - groupMetastatic, levels=design)
res_hlt_vs_met <- test_pseudobulk_fixed(fit, contrast_healthy_metastatic)

stats <- res_met_vs_hlt$logFC
names(stats) <- rownames(res_met_vs_hlt)

stats <- res_hlt_vs_met$logFC
names(stats) <- rownames(res_hlt_vs_met)

#Run GSEA
set.seed(1234)
fgsea_res<- fgsea(pathways = list(pro_CTC_signature = up_CTC),
                  stats = stats)

#AUCELL
expr_matrix <- bulk_expr

cell_ranking <- AUCell_buildRankings(as.matrix(expr_matrix), plotStats = FALSE)
gene_sets <- list(pro_CTC_signature= up_CTC)

cell_AUC <- AUCell_calcAUC(gene_sets, cell_ranking)

auc_scores <- as.data.frame(t(getAUC(cell_AUC)))
auc_scores$sample <- rownames(auc_scores)
auc_scores$condition <- ifelse(grepl("MV", auc_scores$sample), "Metastatic", "Healthy")

library(ggplot2)
library(ggpubr)

aucell_plot<- ggplot(auc_scores, aes(x = condition, y = pro_CTC_signature, fill = condition)) +
  geom_boxplot() +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(
    title = "AUC Scores for pro-CTC Signature lung vs CTC",
    x = "Condition",
    y = "AUC Score"
  ) +
  theme_minimal()
ggsave(filename = "/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/CTCs/pseudobulk/aucell_up_lung_vs_ctc.png", plot=aucell_plot,width= 11, height= 9,dpi=300)

#ssGSEA
up_CTC_filtered <- intersect(up_CTC, rownames(bulk_expr))
gene_sets <- list(pro_CTC_signature = up_CTC_filtered)
bulk_matrix <- as.matrix(bulk_expr)

gsvaPar <- ssgseaParam( bulk_matrix, gene_sets)
ssgsea_scores <- gsva(gsvaPar, verbose=FALSE)

ssgsea_df <- as.data.frame(t(ssgsea_scores))
ssgsea_df$sample <- rownames(ssgsea_df)
ssgsea_df$condition <- ifelse(grepl("MV", ssgsea_df$sample), "Metastatic", "Healthy")
wilcox.test(pro_CTC_signature ~ condition, data = ssgsea_df)

# Viz
ssgsea_plot<- ggplot(ssgsea_df, aes(x = condition, y = pro_CTC_signature, fill = condition)) +
  geom_boxplot() +
  labs(
    title = "ssGSEA scores for pro-CTC signature lung vs CTC",
    x = "Condition",
    y = "ssGSEA Enrichment Score"
  ) +
  stat_compare_means(method = "wilcox.test", label.y = max(ssgsea_df$pro_CTC_signature) * 1.1) +
  theme_minimal()
ggsave(filename = "/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/CTCs/pseudobulk/ssgsea_up_lung_vs_ctc.png", plot=ssgsea_plot,width= 11, height= 9,dpi=300)

# check expression levels of up_CTC genes in bulk
up_CTC_filtered <- intersect(up_CTC, rownames(bulk_expr))
signature_expr <- bulk_expr[up_CTC_filtered, ]

library(pheatmap)

# Group info for annotation
annotation_col <- data.frame(Condition = factor(ifelse(grepl("MV", colnames(signature_expr)), "Metastatic", "Healthy")))
rownames(annotation_col) <- colnames(signature_expr)

# Plot heatmap
heatmap_plot <- pheatmap::pheatmap(signature_expr,
                                   cluster_cols = TRUE,
                                   scale = "row",
                                   annotation_col = annotation_col,
                                   show_rownames = FALSE,
                                   show_colnames = TRUE,
                                   main = "Expression of up-CTC Signature in Bulk Samples")

ggsave(filename = "/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/CTCs/pseudobulk/heatmap_up_ctc.png", plot=heatmap_plot,width= 11, height= 9,dpi=300)



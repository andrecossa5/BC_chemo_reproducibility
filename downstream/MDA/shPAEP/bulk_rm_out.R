#shPAEP bulk analysis with EdgeR

library(tximport)
library(DESeq2)
library(readr)
library(tidyr)
library(dplyr)
library(sva)
library(edgeR)
library(tibble)
library(ggplot2)
library(ggfortify)
library(limma)
library(openxlsx)
library(readxl)
library(ggrepel)

#removed shPAEP2_1 and shPAEP2_2
sample<- factor(c("PT_shSCR_1","PT_shSCR_2","PT_shSCR_3",  
                  "PT_shSCR_4","PT_shSCR_5","PT_shPAEP1_1",  
                  "PT_shPAEP1_2","PT_shPAEP1_3","PT_shPAEP1_4",
                  "PT_shPAEP1_5",
                  "PT_shPAEP2_3","PT_shPAEP2_4")) 

#3 condition SCR, PAEP1 & PAEP2
condition_2 <- factor(ifelse(grepl("shPAEP1", sample), "PT_shPAEP1",
                           ifelse(grepl("shPAEP2", sample), "PT_shPAEP2",
                                  "PT_shSCR")),
                    levels  = c("PT_shSCR", "PT_shPAEP1", "PT_shPAEP2"))

#2 condition SCR & PAEP
condition <- factor(ifelse(grepl("shPAEP1", sample), "PT_shPAEP",
                             ifelse(grepl("shPAEP2", sample), "PT_shPAEP",
                                    "PT_shSCR")),
                      levels = c("PT_shSCR", "PT_shPAEP"))

cell_cycle<- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4",
               "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL",     
               "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2",   
               "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2",
               "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
               "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
               "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1",
               "E2F8", "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
               "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
               "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
               "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",  
               "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
               "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
               "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
               "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

sampletable <- data.frame(
  sample = sample,
  condition = condition
)

files<- file.path(sample, "quant.sf")
names(files)<- sample
tx2gene<- read.delim("tx2gene.tsv", header=FALSE, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("transcript_id", "gene_id", "gene_name")
txi<-tximport(files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

#Create dds with raw counts
dds <- DESeqDataSetFromTximport(txi, sampletable, ~condition)
counts<- dds@assays@data@listData[["counts"]]

#rename ENSG with gene ID  
gene_id_to_name <- tx2gene %>% distinct(gene_id, gene_name) %>% column_to_rownames("gene_id")
rownames(counts) <- gene_id_to_name[rownames(counts), "gene_name"]
counts <- counts[!is.na(rownames(counts)), ]

#Create DGElist
y <- DGEList(counts = counts, group = condition)



keep <- filterByExpr(y)
y<- y[keep, ,keep.lib.sizes=FALSE]
y<- calcNormFactors(y)
design<-model.matrix(~ 0 + condition)
y<- estimateDisp(y,design=design)

### PCA ###
cv_function <- function(x) {
  if (abs(mean(x)) <= 1e-28) return(0)  
  return(sd(x) / mean(x))
}
st_var <- function(x) {
  (x - mean(x)) / sd(x)
}


logCPM <- t(cpm(y, log=TRUE, prior.count = 2))
logCPM_t<- t(logCPM)

#library_size values & cell_cycle_mean values
num_reads<- colSums(y$counts)
cell_cycle_reads <- logCPM_t[rownames(logCPM_t) %in% cell_cycle, ]
mean_cell_cycle_reads <- colMeans(cell_cycle_reads, na.rm=TRUE)

### no correction ###
cv_values_cpm <- apply(logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

#check mean and var # 
mean_st_var_values_cpm <- apply(st_var_values_cpm, 2, mean)
var_st_var_values_cpm <- apply(st_var_values_cpm, 2, var)

pca<-prcomp(st_var_values_cpm)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)
pca_var <- summary(pca)$importance[2, ] 
pc1_var <- pca_var[1]  
pc2_var <- pca_var[2]
condition<- condition_2
pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 1000") +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc1_var * 100, 2), "% variance)"),  
       y = paste0("PC2 (", round(pc2_var * 100, 2), "% variance)")) +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
         label = paste("PC1 vs Cell Cycle Correlation: ",
                       round(cor_pc1_cellcycle, 2)),
         hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/res_final/pca_1000_rm_out.png", plot=pca_plot,width= 9, height= 6,dpi=300)

### regression num_reads
remove_num_reads_effect <- function(gene_expr) {
  model <- lm(gene_expr ~ num_reads)
  residuals(model)  
}

corrected_logCPM_num <- apply(logCPM, 2, remove_num_reads_effect)
cv_values_cpm <- apply(corrected_logCPM_num, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

pca<-prcomp(st_var_values_cpm, scale. = TRUE)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)
pca_var <- summary(pca)$importance[2, ] 
pc1_var <- pca_var[1]  
pc2_var <- pca_var[2]
pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 1000") +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc1_var * 100, 2), "% variance)"),  
       y = paste0("PC2 (", round(pc2_var * 100, 2), "% variance)")) +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")
ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/res_final/pca_1000_num_reads_rm_out.png", plot=pca_plot,width= 9, height= 6,dpi=300)

### regression cell_cycle ###
remove_cell_cycle_effect <- function(gene_expr) {
  model <- lm(gene_expr ~ mean_cell_cycle_reads)
  residuals(model)  
}

corrected_logCPM <- apply(logCPM, 2, remove_cell_cycle_effect)
cv_values_cpm <- apply(corrected_logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

pca<-prcomp(st_var_values_cpm, scale. = TRUE)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)
pca_var <- summary(pca)$importance[2, ] 
pc1_var <- pca_var[1]  
pc2_var <- pca_var[2]
pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition_2)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 1000") +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc1_var * 100, 2), "% variance)"),  
       y = paste0("PC2 (", round(pc2_var * 100, 2), "% variance)")) +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/res_final/pca_1000_cell_cycle_rm_out.png", plot=pca_plot,width= 9, height= 6,dpi=300)

### combat ###
num_reads_scaled <- scale(num_reads)
mod <- model.matrix(~ num_reads_scaled + mean_cell_cycle_reads)
batch <- factor(c(rep("PT_shSCR", 5), rep("PT_shPAEP1", 5), rep("PT_shPAEP2", 2)))
corrected_logCPM <- ComBat(dat = as.matrix(logCPM_t), batch = batch , mod = mod)

corrected_logCPM_t<-t(corrected_logCPM)
cv_values_cpm <- apply(corrected_logCPM_t, 2, cv_function)

cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)


pca<-prcomp(st_var_values_cpm, scale. = TRUE) 
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)
pca_var <- summary(pca)$importance[2, ] 
pc1_var <- pca_var[1]  
pc2_var <- pca_var[2]
pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition_2)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 1000") +
  theme_bw() +
  labs(x = paste0("PC1 (", round(pc1_var * 100, 2), "% variance)"),  
       y = paste0("PC2 (", round(pc2_var * 100, 2), "% variance)")) +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/res_final/pca_1000_combat_rm_out.png", plot=pca_plot,width= 9, height= 6,dpi=300)

### Deseq2 ###
# smallestGroupSize <- 3
# keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
# 
# gene_id_to_name <- tx2gene %>% distinct(gene_id, gene_name) %>% column_to_rownames("gene_id")
# gene_ids <- rownames(dds)
# gene_names <- gene_id_to_name[gene_ids, "gene_name"]
# valid_genes <- !is.na(gene_names)
# dds <- dds[valid_genes, ]
# rownames(dds) <- gene_names[valid_genes]
# 
# dds_filtered <- dds[keep,]
# dds_filtered$condition <- factor(dds_filtered$condition, levels = c(levels(dds_filtered$condition), "PT_shPAEP"))
# dds_filtered$condition[dds_filtered$condition %in% c("PT_shPAEP1", "PT_shPAEP2")] <- "PT_shPAEP"
# dds_filtered$condition <- droplevels(dds_filtered$condition)
# dds_final <- DESeq(dds_filtered)
# 
# 
# 
# 
# rld <- rlog(dds, blind = FALSE)
# vsd <- vst(dds, blind = FALSE)
# pca_plot<- plotPCA(rld, intgroup=c("condition"))
# pca_plot_vsd<- plotPCA(vsd, intgroup=c("condition"))
# 
# pca_plot$data
# pca_plot_vsd$data
# 
# ### deseq2 analysis ###
# dds_final$condition
# res <- results(dds_final, contrast=c("condition","PT_shSCR","PT_shPAEP"))
# summary(res)
# resOrdered <- res[order(res$pvalue),]
# res_filtered<-subset(resOrdered, padj < 0.1)
# 
# 

### DE edgeR ###
fit <- glmQLFit(y, design)

SCRvsPAEP <- makeContrasts(conditionPT_shPAEP - conditionPT_shSCR,
                           levels = design)

qlf<-glmQLFTest(fit, contrast=SCRvsPAEP)
t<-topTags(qlf, n=Inf)
tt<- t$table

ttt<-tt[tt$FDR<=0.1,]
ttt<- ttt %>% arrange(desc(logFC))
write.xlsx(genes_deg,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_paep_vs_scr_rm_out_thrb_regulon.xlsx",rowNames = TRUE)

#genes regulon THRB 
genes<- c('TP63', 'EXD2', 'USP44', 'CHI3L2', 'RAB33A', 'THRB', 'UCP3', 
           'HMGCLL1', 'ATP2A3', 'ANKS1A', 'GBX2', 'SCN2A', 'PAEP', 'HR', 
           'ARHGAP31', 'FHDC1', 'ANGPT1', 'NOD1')
logCPM <- t(cpm(y, log=TRUE, prior.count = 2))
logCPM_t<- t(logCPM)
write.csv(logCPM_t,'/Users/ieo7295/Desktop/BC_chemo_reproducibility/data/MDA/bulk_expr_matr.csv')

#DE regulon genes in shPAEP vs shSCR
fit <- glmQLFit(y, design)

SCRvsPAEP <- makeContrasts(conditionPT_shPAEP - conditionPT_shSCR,
                           levels = design)

qlf<-glmQLFTest(fit, contrast=SCRvsPAEP)
t<-topTags(qlf, n=Inf)
tt<- t$table

ttt<-tt[tt$FDR<=0.1,]
ttt<- ttt %>% arrange(desc(logFC))

genes_deg<- tt[rownames(tt) %in% genes, ]
sample_scores <- colMeans(genes_deg)

sample_info <- data.frame(
  Sample = names(sample_scores),
  ActivityScore = sample_scores,
  Condition = ifelse(grepl("shSCR", names(sample_scores)), "shSCR", "shPAEP")
)
write.csv(sample_info,'/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/MDA/thrb_genes_expression_bulk.csv')

# Plot: Boxplot + Wilcoxon
library(ggpubr)
ggboxplot(sample_info, x = "Condition", y = "ActivityScore", add = "jitter") +
  stat_compare_means(method = "wilcox.test") +
  theme_minimal()



#using pca to give THRB regulon activity score
pca_result <- prcomp(t(heat_data))
regulon_activity_scores_pca <- pca_result$x[, 1]  
regulon_activity_scores_pca_normalized <- scale(regulon_activity_scores_pca)
activity_data <- data.frame(ActivityScore = regulon_activity_scores_pca_normalized, Condition = condition_2)
ggplot(activity_data, aes(x = condition_2, y = ActivityScore)) +
  geom_boxplot() +
  labs(x = "Condition", y = "Regulon Activity Score (Normalized)") +
  theme_minimal()

#Regulon enrichment




library(ggpubr)

my_comparisons <- list(
  c("shSCR", "shPAEP1"),
  c("shSCR", "shPAEP2"),
  c("shPAEP1", "shPAEP2")
)
activity_data$Condition <- factor(activity_data$Condition, levels = c("shSCR", "shPAEP1", "shPAEP2"))
write.csv(,"/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/MDA/enrichment_thrb_sh.csv")



#PAEPvsSCR
PAEPvsSCR<- makeContrasts(conditionPT_shPAEP - conditionPT_shSCR,
                          levels=design)
qlf<-glmQLFTest(fit, contrast=PAEPvsSCR)
k<-topTags(qlf, n=Inf)
kk<-k$table

kkk<-kk[kk$FDR<=0.1,]
kkk<- kkk %>% arrange(desc(logFC))
write.xlsx(kk,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_scr_vs_paep_rm_out.xlsx",rowNames = TRUE)

### PAEP1 vs PAEP2 ### 
y_3 <- DGEList(counts = counts, group = condition_2)


keep <- filterByExpr(y_3)
y_3<- y_3[keep, ,keep.lib.sizes=FALSE]
y_3<- calcNormFactors(y_3)
design_3<-model.matrix(~ 0 + condition_2)
y_3<- estimateDisp(y,design=design_3)

fit <- glmQLFit(y_3, design_3)

PAEP1vsPAEP2 <- makeContrasts(condition_2PT_shPAEP1 - condition_2PT_shPAEP2,
                              levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsPAEP2)
x_3<-topTags(qlf, n=Inf)
xx_3<- x_3$table

xxx_3<-xx_3[xx_3$FDR<=0.1,]
xxx_3<-xxx_3 %>% arrange(desc(logFC))
write.xlsx(kkk,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_paep1_vs_paep2_rm_out.xlsx",rowNames = TRUE)

### PAEP1 vs SCR ###
PAEP1vsSCR <- makeContrasts(condition_2PT_shPAEP1 - condition_2PT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsSCR)
t_3<-topTags(qlf, n=Inf)
tt_3<- t_3$table

ttt_3<-tt_3[tt_3$FDR<=0.1,]
ttt_3<-ttt_3 %>% arrange(desc(logFC))
write.xlsx(ttt_3,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_paep1_vs_scr_rm_out.xlsx",rowNames = TRUE)

### PAEP2vsSCR ###
PAEP2vsSCR <- makeContrasts(condition_2PT_shPAEP2 - condition_2PT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP2vsSCR)
n_3<-topTags(qlf, n=Inf)
nn_3<- n_3$table

nnn_3<-nn_3[nn_3$FDR<=0.1,]
nnn_3<-nnn_3 %>% arrange(desc(logFC))
write.xlsx(nnn_3,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_paep2_vs_scr_rm_out.xlsx",rowNames = TRUE)


### GRN all contrasts###
grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr <- tt[rownames(tt) %in% grn_corr$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.05) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor <- tt[rownames(tt) %in% grn_chemor$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.05) %>% 
  arrange(desc(logFC)) 
write.xlsx(grn_ensg_corr,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_corr_rm_out.xlsx",rowNames = TRUE)
write.xlsx(grn_ensg_chemor,"/Users/ieo7295/Desktop/BC_sh/results/res_final/Degs_chemor_rm_out.xlsx",rowNames = TRUE)

### heatmap GRN ###
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(grn_ensg_chemor), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(heat_data,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_2, row.names = colnames(heat_data)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette)

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/res_final/heatmap_grn_chemor.png", plot=heatmap_plot,width= 9, height= 6,dpi=300)

### volcano_plot ####
volcano_plot <- function(ttt, title) {
  
  ttt$logFC <- as.numeric(as.character(ttt$logFC))
  ttt$threshold <- ttt$FDR < 0.1 & abs(ttt$logFC) > 1
  
  ttt$Gene <- rownames(ttt)
  
  ggplot(ttt, aes(x = logFC, y = -log10(PValue), color = threshold)) + 
    geom_point(alpha = 0.8) + 
    scale_color_manual(values = c("black", "red")) + 
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),  
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5), 
      panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25)  
    ) + 
    ggtitle(title) + 
    xlab("Log Fold Change") + 
    ylab("-Log10 P-value") + 
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) + 
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_text_repel(data = subset(ttt,threshold), aes(label = Gene), 
                    size = 3, max.overlaps = 25, 
                    box.padding = 0.9, point.padding = 0.5,
                    force=2,
                    nudge_y= 0.5,
                    vjust=1,
                    color="black")
}

volcano_plot <- volcano_plot(tt,title ='Differential Gene Expression')
ggsave(plot= volcano_plot, "/Users/ieo7295/Desktop/BC_sh/results/res_final/volcano_plot_rm_out.png", dpi=300)

# Heatmap plot
top_genes<- tail(ttt,50)
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(top_genes), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot <- pheatmap::pheatmap(heat_data,
                         cluster_cols = TRUE,
                         scale = "row",
                         annotation_col = data.frame(condition = condition_2, row.names = colnames(logCPM_t)),
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         color = color_palette)

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/res_final/heatmap_paep_vs_scr_rm_out_bottom50.png", plot=heatmap_plot,width= 11, height= 9,dpi=300)

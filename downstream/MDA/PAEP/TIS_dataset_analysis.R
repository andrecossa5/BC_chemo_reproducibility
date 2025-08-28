#!/usr/bin/env Rscript

# script to analyze PAEP expression in TIS dataset (Bajtai et al.,2025)

library(edgeR)
library(limma)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(biomaRt)

# ====== INPUTS ======
count_file <- "GSE287953_count_matrix_bulk.txt"  # path to your count matrix
target_gene <- "PAEP"
readLines(count_file, n = 20)
# ====== READ COUNTS ======
counts <- read.table(count_file, header=TRUE, sep=" ", row.names=1, check.names=FALSE)


if (!requireNamespace("biomaRt", quietly=TRUE)) install.packages("biomaRt")

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_map <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                  filters="ensembl_gene_id",
                  values=sub("\\..*$","", rownames(counts)),
                  mart=mart)

# Merge mapping into counts
counts$ensembl_gene_id <- sub("\\..*$","", rownames(counts))
counts_annot <- merge(gene_map, counts, by="ensembl_gene_id")
library(org.Hs.eg.db)
AnnotationDbi::select(org.Hs.eg.db, 
                      keys = "PAEP", 
                      keytype = "SYMBOL", 
                      columns = "ENSEMBL")
# Filter for PAEP
paep_counts <- counts_annot[counts_annot$hgnc_symbol == "PAEP", ] #Hugo name coding

paep_ensg <- "ENSG00000122133"
paep_row <- counts[ sub("\\..*", "", rownames(counts)) == paep_ensg, ] #ensembl 
rownames(counts) <- sub("\\..*$", "", rownames(counts))
# ====== EXTRACT SAMPLE METADATA FROM COLUMN NAMES ======

meta <- data.frame(sample=colnames(counts))
meta$group <- ifelse(grepl("ctr", meta$sample, ignore.case=TRUE), "CTR",
                     ifelse(grepl("tis", meta$sample, ignore.case=TRUE), "TIS",
                            ifelse(grepl("repop", meta$sample, ignore.case=TRUE), "REPOP", NA)))

meta$group <- factor(meta$group, levels=c("CTR","TIS","REPOP"))

meta <- meta %>%
  mutate(
    cell_line = sub("-.*","", sample),
    treatment = group
  )


# ====== FILTER TO MDA231 ======
meta <- meta %>% filter(grepl("231", sample))
counts <- counts[, meta$sample]

# ====== SET FACTORS ======
meta$treatment <- factor(meta$treatment, levels=c("CTR","TIS","REPOP"))

# ====== LIMMA-VOOM ANALYSIS ======
dge <- DGEList(counts=counts)
keep <- filterByExpr(dge, group=meta$group)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~0 + meta$group, data=meta)
colnames(design) <- levels(meta$group)

v <- voom(dge, design, plot=FALSE)
contr <- makeContrasts(
  TIS_vs_CTR = TIS - CTR,
  REPOP_vs_TIS = REPOP - TIS,
  levels=design
)
fit <- lmFit(v, design)
fit <- eBayes(contrasts.fit(fit, contr))

# ====== GET PAEP RESULTS ======
target_gene <- "ENSG00000122133"   

# ===== Extract PAEP from contrasts =====
paep_res <- lapply(colnames(contr), function(cn) {
  tt <- topTable(fit, coef=cn, number=Inf, sort.by="none")
  
  if ("Gene_ID" %in% colnames(tt)) {
    tt[tt$Gene_ID == target_gene, , drop=FALSE]
  } else {
    tt[rownames(tt) == target_gene, , drop=FALSE]
  }
})
names(paep_res) <- colnames(contr)
paep_df <- bind_rows(paep_res, .id="contrast")
print(paep_df)

# ====== BOX PLOT ======
paep_expr <- data.frame(logCPM = v$E[target_gene, ], treatment=meta$group)
p<- ggplot(paep_expr, aes(x=treatment, y=logCPM)) +
  geom_boxplot() + geom_jitter(width=0.1) +
  labs(title=paste("PAEP expression in MDA231"), y="log2 CPM") +
  theme_bw()

ggsave("/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/MDA/TIS/PAEP_expression_MDA231.png", plot=p, width=6, height=4, dpi=300)



# ===== INPUTS =====
target_gene <- "ENSG00000122133"  
contrasts <- colnames(contr)     

# ===== FUNCTION: Get DEG table =====
get_deg_table <- function(fit, coef, mart) {
  tt <- topTable(fit, coef=coef, number=Inf, sort.by="P")
  
  # Add Ensembl IDs as column if only rownames
  if (!("ensembl_gene_id" %in% colnames(tt))) {
    tt <- tt %>% tibble::rownames_to_column("ensembl_gene_id")
  }
  
  # Map Ensembl → HGNC
  gene_map <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
                    filters="ensembl_gene_id",
                    values=sub("\\..*$","", tt$ensembl_gene_id),
                    mart=mart)
  
  tt <- tt %>%
    left_join(gene_map, by="ensembl_gene_id") %>%
    relocate(hgnc_symbol, .before=logFC)
  
  return(tt)
}
if (!requireNamespace("openxlsx", quietly=TRUE)) install.packages("openxlsx")
library(openxlsx)
# ===== Prepare mart =====
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
outdir<- "/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/MDA/TIS"
# ===== Loop through contrasts =====
# Define PAEP ENSG
paep_id <- "ENSG00000122133"

for (cn in contrasts) {
  cat("Processing contrast:", cn, "\n")
  
  # Get annotated DEG list
  tt <- get_deg_table(fit, cn, mart)
  tt <- as.data.frame(tt) %>% tibble::as_tibble()
  
  # Save CSV
  
  out_xlsx <- file.path("/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/MDA/TIS",
                        paste0("DEG_", cn, ".xlsx"))
  
  openxlsx::write.xlsx(tt, out_xlsx, row.names=FALSE)
  cat("Saved:", out_xlsx, "\n")
  
  # Mark PAEP
  tt <- tt %>% mutate(
    is_paep = (ensembl_gene_id == paep_id),
    sig = adj.P.Val < 0.05 & abs(logFC) > 1
  )
  
  # Top 15 genes by adjusted p-value
  top15 <- tt %>% arrange(adj.P.Val) %>% dplyr::slice(1:15)
  
  p <- ggplot(tt, aes(x=logFC, y=-log10(P.Value))) +
    geom_point(aes(color = case_when(
      is_paep ~ "PAEP",
      sig ~ "Significant",
      TRUE ~ "NS"
    )), alpha=0.7) +
    scale_color_manual(
      values=c("NS"="grey70", "Significant"="red", "PAEP"="blue"),
      name=NULL   # removes "case_when" legend title
    ) +
    geom_vline(xintercept=c(-1,1), linetype="dashed") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    geom_text_repel(data=top15, aes(label=hgnc_symbol),
                    size=3, max.overlaps=Inf) +
    geom_text_repel(data=tt %>% filter(is_paep),
                    aes(label=hgnc_symbol),
                    color="blue", size=3, max.overlaps=Inf) +
    labs(title=paste("Volcano plot:", cn),
         x="log2 Fold Change", y="-log10 P-value") +
    theme_bw()
  
  out_png <- file.path(outdir, paste0("Volcano_", cn, "_PAEP.png"))
  ggsave(out_png, p, width=6, height=5, dpi=300)
  cat("Saved:", out_png, "\n")
}

summary(v$E["ENSG00000122133", ])
which(rownames(v$E) == "ENSG00000122133") %in% which(keep)
cpm_vals <- cpm(dge, log=TRUE, prior.count=1)["ENSG00000122133", ]
"ENSG00000122133" %in% rownames(dge)
"ENSG00000122133" %in% rownames(v$E)



paep_expr <- data.frame(
  logCPM = v$E[target_gene, ],
  treatment = meta$group
)

p <- ggplot(paep_expr, aes(x=treatment, y=logCPM)) +
  geom_boxplot() +
  geom_jitter(width=0.1) +
  labs(title="PAEP expression in MDA231", y="log2 CPM") +
  theme_bw() +
  stat_compare_means(
    comparisons = list(c("CTR","TIS"), c("TIS","REPOP"), c("CTR","REPOP")),
    method = "wilcox.test"
  )



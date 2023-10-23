# DE with edgeR
library(tidyverse)
library(SingleCellExperiment)
library(zellkonverter)
library(edgeR)


##


# Utils
aggregate_data <- function(sce, assay='raw', group_by=NULL, cont_covariates=c('nUMIs', 'mito_perc'), 
                           cat_covariates=c('seq_run')
  ) {
  
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


##


fit_pseudobulk_fixed <- function(df_pseudo, sce, test_column=NULL, design=NULL) {
  
  # Load DGEList
  y <- DGEList(
    df_pseudo[, colnames(df_pseudo) %in% row.names(sce)] %>% t(), 
    samples=df_pseudo[, !colnames(df_pseudo) %in% row.names(sce)],
    group=test_column
  )
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  # Estimate dispersion and fit coefficient glm
  y <- estimateDisp(y, design=design)
  fit <- glmQLFit(y, design)
  
  return(fit)

}


##


fit_pseudobulk_random <- function(df_pseudo, sce, test_column=NULL, block_column=NULL, design=NULL) {
  
  # Load DGEList
  y <- DGEList(
    df_pseudo[, colnames(df_pseudo) %in% row.names(sce)] %>% t(), 
    samples=df_pseudo[, !colnames(df_pseudo) %in% row.names(sce)],
    group=test_column
  )
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  # Estimate dispersion and fit coefficient with random effects
  fit <- voomLmFit(y, design=design, block=block_column, sample.weights=TRUE)
  
  return(fit)
  
}


##


test_pseudobulk_fixed <- function(fit, contrast=NULL) {
  
  # Test contrasts
  qlf <- glmQLFTest(fit, contrast=contrast)
  
  tt <- topTags(qlf, n=Inf)
  tt <- tt$table
  tt <- tt %>% arrange(desc(logFC))
  
  return(tt)
  
}


##



test_pseudobulk_random <- function(fit, contrast=NULL) {
  
  # Test contrasts
  con.fit <- contrasts.fit(fit, contrasts=contrast) 
  con.fit <- eBayes(con.fit)
  
  tt <- topTable(con.fit, n=Inf)
  tt <- tt %>% arrange(desc(logFC))
  
  return(tt)
  
}


##


# Paths
path_main <- '/Users/IEO5505/Desktop/BC_chemo_reproducibility/'
path_data <- paste0(path_main, '/data/MDA')
path_results <- paste0(path_main, '/results/MDA/pseudobulk')

# Load adata as sce
sce <- readH5AD(paste0(path_data, '/clustered.h5ad'))


##


# Pseudo-bulk

# Filter GBC-sample combinations with less than 10 cells
min_cells <- 10
colData(sce)$group <- paste0(colData(sce)$GBC, colData(sce)$sample)
filtered_sce <- sce[, colData(sce)$group %in% names(which(table(colData(sce)$group) >= min_cells)) ]
colData(filtered_sce)$group <- as.factor(colData(filtered_sce)$group)

# Aggregate into a pseudobulk df. For continuous covariates, the median is used, for GE, the sum.
df_pseudo <- aggregate_data(filtered_sce, assay='raw', group_by='group', 
                            cont_covariates=c('nUMIs', 'mito_perc', 'G1.S', 'G2.M'),
                            cat_covariates=c('seq_run', 'condition', 'dataset', 'origin'))

# Psudobulk DE with edgeR
df_pseudo[,1:10] %>% head()
df_pseudo %>% dim()

# Here we go

# Model specification ~ 0 + branch_origin + confounders
batch <- as.factor(as.character(df_pseudo$seq_run))
nUMIs <- df_pseudo$nUMIs
mito_perc <- df_pseudo$mito_perc
G1.S <- df_pseudo$G1.S
G2.M <- df_pseudo$G2.M
block_column <- df_pseudo$dataset

test_column <- as.factor(paste0(df_pseudo$condition, '_', df_pseudo$origin))
design <- model.matrix(~ 0 + test_column + nUMIs + mito_perc)
design %>% head()
confounders <- '_nUMIs_mito_perc'

# Create contrasts
PTs_vs_lungs <- makeContrasts(
  (test_columnAC_AC_PT+test_columnAC_NT_PT+test_columnNT_AC_PT+test_columnNT_NT_PT)/4-(test_columnAC_AC_lung+test_columnAC_NT_lung+test_columnNT_AC_lung+test_columnNT_NT_lung)/4,
  levels=design
)
single_treated_vs_untreated_lungs <- makeContrasts(test_columnAC_NT_lung-test_columnNT_NT_lung, levels=design)
double_treated_vs_untreated_lungs <- makeContrasts(test_columnAC_AC_lung-test_columnNT_NT_lung, levels=design)
double_treated_vs_single_treated_lungs <- makeContrasts(test_columnAC_AC_lung-test_columnAC_NT_lung, levels=design)
treated_vs_untreated_PTs <- makeContrasts((test_columnAC_AC_PT+test_columnAC_NT_PT)/2-(test_columnNT_NT_PT+test_columnNT_AC_PT)/2, levels=design)
myContrasts <- list(
  PTs_vs_lungs, single_treated_vs_untreated_lungs, double_treated_vs_untreated_lungs, 
  double_treated_vs_single_treated_lungs, treated_vs_untreated_PTs
)

# Fit fixed and random models
fit_fixed <- fit_pseudobulk_fixed(df_pseudo, filtered_sce, test_column, design)
fit_random <- fit_pseudobulk_random(df_pseudo, filtered_sce, test_column, block_column, design)

# Test contrasts
names_contrasts <- c('PTs_vs_lungs', 'single_treated_vs_untreated_lungs', 'double_treated_vs_untreated_lungs', 
                     'double_treated_vs_single_treated_lungs', 'treated_vs_untreated_PTs')
res_fixed <- lapply(myContrasts, function(x) { test_pseudobulk_fixed(fit_fixed, x) } )
res_fixed <- setNames(res_fixed, names_contrasts)
res_fixed <- lapply(myContrasts, function(x) { test_pseudobulk_random(fit, x) } )


# Save results fixed
for (x in names(res_fixed)) {
  df <- res_fixed[[x]] %>% mutate(contrast=x) 
  write.csv(df, paste0(path_results, '/fixed/results_', x, confounders, '.csv'))
}

# Save results random
for (x in names(res_random)) {
  df <- res_random[[x]] %>% mutate(contrast=x) 
  write.csv(df, paste0(path_results, '/random/results_', x, confounders, 'dataset.csv'))
}


##




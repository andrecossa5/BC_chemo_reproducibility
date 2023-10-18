# Prep reference from Pai et al., 2021

library(tidyverse)
library(Seurat)
library(data.table)


##

# Path data
path_data <- '/Users/IEO5505/Desktop/BC_chemo_reproducibility/data/Pai_reference'

# Read files
rds_files <- list.files(path_data, pattern="\\.rds$", full.names=TRUE)
data_list <- lapply(rds_files, readRDS)
names(data_list) <- c('TNBC', 'TNBCSub', 'TNBCTC', 'TNBCTum')


data_list$TNBC@meta.data %>% colnames()
data_list$TNBC@meta.data$group %>% unique()
data_list$TNBC@meta.data$seurat_clusters %>% unique()


glimpse(data_list$TNBC)


all((data_list$TNBCTum %>% colnames()) %in% (data_list$TNBC %>% colnames()))

data_list$TNBC@meta.data$seurat_clusters %>% unique()
data_list$TNBCTum@meta.data$seurat_clusters %>% unique()

M <- data_list$TNBC@assays$RNA@counts %>% as.matrix() %>% t()


DimPlot(data_list$TNBCTum)




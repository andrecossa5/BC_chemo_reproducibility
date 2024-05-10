
### Run Lin, Poi or Nb regression ###

# Notebook with trajectory analysis at: https://github.com/sankaranlab/redeem_reproducibility/blob/master/Note-6%20HSC%20clonal%20output%20analysis.ipynb
# Script to run Regression at: https://github.com/sankaranlab/redeem_reproducibility/blob/master/API/Run_Lin_regression.R

library(dplyr)
library(tibble)
library(openxlsx)
library(ggrepel)

SEED <- 4321
set.seed(SEED)

path_main <- "/Users/ieo6983/Desktop/breast_alberto/"
path_utils <- fs::path(path_main, "BC_chemo_reproducibility/downstream/MDA/trajectory/")
path_input_data <- fs::path(path_main, "data/MDA/")
path_results <- fs::path(path_main, "results/MDA")

if(!dir.exists(path_results)){
  dir.create(path_results, recursive = T)
}

source(fs::path(path_utils,"utils_trajectory.R"))


##


args = commandArgs(trailingOnly=TRUE)

if(!is_empty(args)){
  mode <- args[1] ## "lm" "poi" "nb"
  core <- as.numeric(args[2]) # Number of cores, I.e. 8
  LinOut.df_rds <- args[3] # path to input data.frame for regression (or file name, if located in 'current')
  name <- args[4] # name for output file with DEGs resulting from regression
} else {
  mode <- "poi"
  core <- as.numeric(8) 
  LinOut.df_rds <- fs::path(path_input_data, "agg_for_poisson.RDS") 
  name <- "reg_out" 
}

## Read in the LinOut.gene.df file
print("Read in")
print(LinOut.df_rds)
LinOut.df<-readRDS(LinOut.df_rds)

boxplot(log(LinOut.df[, c(14:44)]))

# Filter out genes with 0-expression in more than 10% of samples
LinOut.df.filt <- filter_genes(LinOut.df)
boxplot(log(LinOut.df.filt[, c(14:44)]))

## Run regression
if(mode=="lm"){
  print("run linear model")
  LinOut.result <- Run_Lin_regression(LinOut.df.filt, n.cores = core)
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result <- Run_Lin_regression_poi(LinOut.df.filt, n.cores = core, qval = T, tot_UMIs = T)
}else if(mode=="nb"){
  # W: nb does not support vectors of zeros as response variable (genes with 0 expression across all clones)
  print("run negative-binomial model")
  LinOut.result <- Run_Lin_regression_nb(LinOut.df.filt, n.cores = core, qval = T, tot_UMIs = T)
}


if("qs" %in% names(LinOut.result)){
  LinOut.result.df <- convert_results_to_df(LinOut.result)    
}else{
  LinOut.result.df <- merge(LinOut.result$slopes, LinOut.result$ps, by = "row.names")
  colnames(LinOut.result.df) <- c("genes", "slopes", "ps")  
}

par(mfrow = c(1, 1))
# qline = 0.05; pline = 0.001
p1 <- PlotLinRegress_Vocano(LinOut.result.df, slot = "ps", pline = 0.01)
p1

# Save output 
#saveRDS(LinOut.result, fs::path(path_results, paste0(name, ".", mode, ".RDS")))
#LinOut.result.df %>% arrange(., qs) %>% write.xlsx(., fs::path(path_results, paste0(name, ".", mode, ".xlsx")), rowNames=T)


##


## Repeat regression with normalization & scaling ##

# Filter and normalize clone-level gene expression
LinOut.df.norm.scale <- normalize_aggregated_expr(LinOut.df, scale = T)

# Aggregate gene expression 
boxplot(log(LinOut.df.norm.scale[, c(14:44)]))

# Regression with normalized and scaled GE matrix
mode <- "poi" # "lm" "poi"
core <- as.numeric(8) # Number of cores, I.e. 8

## Run regression
if(mode=="lm"){
  print("run linear model")
  LinOut.result.norm.scale <- Run_Lin_regression(LinOut.df.norm.scale, n.cores=core)
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result.norm.scale <- Run_Lin_regression_poi(LinOut.df.norm.scale, n.cores=core, qval = F, tot_UMIs = T)
}


if("qs" %in% names(LinOut.result)){
  LinOut.result.norm.scale.df <- convert_results_to_df(LinOut.result.norm.scale)    
}else{
  LinOut.result.norm.scale.df <- merge(LinOut.result.norm.scale$slopes, LinOut.result.norm.scale$ps, by = "row.names")
  colnames(LinOut.result.norm.scale.df) <- c("genes", "slopes", "ps")  
}

p2 <- PlotLinRegress_Vocano(LinOut.result.norm.scale.df, slot = "ps")
p2


##


## Repeat regression with normalization  ##

# Filter and normalize clone-level gene expression
LinOut.df.norm <- normalize_aggregated_expr(LinOut.df, scale = F)

# Aggregate gene expression 
boxplot(log(LinOut.df.norm[, c(14:44)]))
boxplot(LinOut.df.norm[, c(14:44)])

# Regression with normalized and scaled GE matrix
mode <- "poi" # "lm" "poi"
core <- as.numeric(8) # Number of cores, I.e. 8

## Run regression
if(mode=="lm"){
  print("run linear model")
  LinOut.result.norm <- Run_Lin_regression(LinOut.df.norm,n.cores=core)
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result.norm <- Run_Lin_regression_poi(LinOut.df.norm, n.cores=core, qval = F, tot_UMIs = T)
}

if("qs" %in% names(LinOut.result.norm)){
  LinOut.result.norm.df <- convert_results_to_df(LinOut.result.norm)    
}else{
  LinOut.result.norm.df <- merge(LinOut.result.norm$slopes, LinOut.result.norm$ps, by = "row.names")
  colnames(LinOut.result.norm.df) <- c("genes", "slopes", "ps")  
}

p3 <- PlotLinRegress_Vocano(LinOut.result.norm.df, slot = "ps")
p3


##


## Run regression based on condition
if(mode=="lm"){
  print("No linear model available for condition")
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result.cond <- Run_Lin_regression_poi.for_condition(LinOut.df.filt, n.cores = core, qval = T, tot_UMIs = T)
}else if(mode=="nb"){
  # W: nb does not support vectors of zeros as response variable (genes with 0 expression across all clones)
  print("run negative-binomial model")
  LinOut.result.cond <- Run_Lin_regression_nb.for_condition(LinOut.df.filt, n.cores = core, qval = T, tot_UMIs = T)
}

LinOut.result.cond.df <- list()
for(exp_cond in names(LinOut.result.cond)){
  print(exp_cond)
  LinOut.result.cond.i <- LinOut.result.cond[[exp_cond]]
  if("qs" %in% names(LinOut.result.cond.i)){
    LinOut.result.cond.i.df <- convert_results_to_df(LinOut.result.cond.i)
    LinOut.result.cond.df[[exp_cond]] <- LinOut.result.cond.i.df  
  }else{
    LinOut.result.cond.i.df <- merge(LinOut.result.cond.i$slopes, LinOut.result.cond.i$ps, by = "row.names") 
    colnames(LinOut.result.cond.i.df) <- c("genes", "slopes", "ps")
    LinOut.result.cond.df[[exp_cond]] <- LinOut.result.cond.i.df  
  }
}
  

p4 <- PlotLinRegress_Vocano(LinOut.result.cond.df$NT_NT, slot = "qs")
print(p4)

p5 <- PlotLinRegress_Vocano(LinOut.result.cond.df$NT_AC, slot = "qs")
print(p5)

p6 <- PlotLinRegress_Vocano(LinOut.result.cond.df$AC_AC, slot = "qs")
print(p6)

# Save output 
#saveRDS(LinOut.result.cond, fs::path(path_results, paste0(name, ".cond", ".", mode, ".RDS")))
#lapply(LinOut.result.cond.df, FUN = function(x){arrange(x, qs)}) %>% write.xlsx(., fs::path(path_results, paste0(name, ".cond", ".", mode, ".xlsx")), rowNames=T)


##


## Data exploration

# Gene_expr distribution
par(mfrow = c(3, 4))

for(i in 1:dim(LinOut.df.filt[, 14:25])[2]){
  gene <- colnames(LinOut.df.filt[13+i])
  
  m <- round(mean(LinOut.df.filt[,13+i]),2)
  v <- round(var(LinOut.df.filt[,13+i]),2)
  
  d <- density(LinOut.df.filt[,13+i])
  p <- plot(d$x, d$y, xlab = paste0("gene-expr counts - ",gene), ylab = "density", 
            main = paste0("Mean: ", m, " Var: ", v))
  print(p)
}

for(i in 1:dim(LinOut.df.filt[, 14:25])[2]){
  print(i)
  gene <- names(LinOut.df.filt)[13+i]
  d <- density(log(LinOut.df.filt[, 13+i]))
  p <- plot(d$x, d$y, xlab = paste0("log(gene-expr counts) - ",gene), ylab = "density")
  print(p)
}

clonesize.umi <- LinOut.df.filt[14:dim(LinOut.df.filt)[2]] %>% rowSums()
LinOut.df.filt.norm <- 1e6*LinOut.df.filt[, 14:dim(LinOut.df.filt)[2]] / clonesize.umi

for(i in 1:dim(LinOut.df.filt.norm[1:12])[2]){
  gene <- colnames(LinOut.df.filt.norm)[i]
  print(gene)
  
  m <- round(mean(LinOut.df.filt.norm[, i]), 2)
  v <- round(var(LinOut.df.filt.norm[,i]),2)
  
  d <- density(LinOut.df.filt.norm[,i])  
  p <- plot(d$x, d$y, xlab = paste0("gene-expr counts - ",gene), ylab = "density", 
            main = paste0("Mean: ", m, " Var: ", v))
  print(p)
}
    

# Gene_expr ~ met_potential
par(mfrow = c(3, 4))

x <- LinOut.df$met_potential
for(i in 1:dim(LinOut.df[, 14:25])[2]){
  y <- LinOut.df[, 14+i]
  p <- plot(y~x)
  print(p)
}

x <- LinOut.df$met_potential
for(i in 1:dim(LinOut.df[, 14:25])[2]){
  y <- log(LinOut.df[, 14+i])
  p <- plot(y~x)
  print(p)
}

par(mfrow = c(1,1))

# Significant p-values
print("Percentage of genes with a significant q-value (all clones):")
print(paste0(round((sum(LinOut.result.df$qs <= 0.05) / dim(LinOut.result.df)[1]) * 100,2), "%"))

for(cond in names(LinOut.result.cond.df)){
  print(paste0("Percentage of genes with a significant p-value (", cond, "):"))
  print(paste0(round((sum(LinOut.result.cond.df[[cond]]$qs <= 0.05) / dim(LinOut.result.cond.df[[cond]])[1]) * 100,2), "%"))  
}




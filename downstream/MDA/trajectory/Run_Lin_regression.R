
### Run Lin or Poi regression ###

# Notebook with trajectory analysis at: https://github.com/sankaranlab/redeem_reproducibility/blob/master/Note-6%20HSC%20clonal%20output%20analysis.ipynb
# Script to run Regression at: https://github.com/sankaranlab/redeem_reproducibility/blob/master/API/Run_Lin_regression.R

library(purrr)
library(dplyr)
library(tibble)
library(openxlsx)
library(ggrepel)

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
  mode <- args[1] ## "lm" "poi"
  core <- as.numeric(args[2]) # Number of cores, I.e. 8
  LinOut.df_rds <- args[3] # path to input data.frame for regression (or file name, if located in 'current')
  name <- args[4] # name for output file with DEGs resulting from regression
} else {
  mode <- "poi" # "lm" "poi"
  core <- as.numeric(8) # Number of cores, I.e. 8
  LinOut.df_rds <- fs::path(path_input_data, "agg_for_poisson.RDS") # path to input data.frame for regression (or file name, if located in 'current')
  name <- "reg_poi.out" # name for output file with DEGs resulting from regression
}

## Read in the LinOut.gene.df file
print("Read in")
print(LinOut.df_rds)
LinOut.df<-readRDS(LinOut.df_rds)

## Run regression
if(mode=="lm"){
  print("run linear model")
  LinOut.result<-Run_Lin_regression(LinOut.df,n.cores=core)
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result<-Run_Lin_regression_poi(LinOut.df,n.cores=core)
}

# Save output 
#saveRDS(LinOut.result, fs::path(path_results, paste0(name,".RDS")))
#write.xlsx(LinOut.result, fs::path(path_results, paste0(name,".xlsx")), rowNames=T)


##

LinOut.result.df <- convert_results_to_df(LinOut.result)

# Volcano plot of results
PlotLinRegress_Vocano <- function(LinOut.result.df, slot = "ps", pline = 0.001, qline = 0.05){
  datatoplot <- LinOut.result.df
  datatoplot <- datatoplot %>% mutate(score = -log(qs)*abs(slopes))
  
  #Label <- subset(datatoplot, qs < 0.2) %>% .[order(.$score,decreasing=T),] %>% .[1:80,] %>% .[, "genes"]
  #slope.L <- subset(datatoplot, qs < 0.2)$slopes %>% min
  #slope.R <- subset(datatoplot,qs < 0.2)$slopes %>% max
  Label <- datatoplot[order(datatoplot$score,decreasing=T),] %>% .[1:80,] %>% .[, "genes"]
  slope.L <- datatoplot$slopes %>% min
  slope.R <- datatoplot$slopes %>% max
  datatoplot$Label <- ifelse(datatoplot$genes %in% Label, datatoplot$genes, "")
  Name = "volcano"
  
  if(slot=="qs"){    
    p <- ggplot(datatoplot)+aes(slopes, -log10(qs), label = Label)+
      geom_point()+
      geom_hline(yintercept=-log10(qline),linetype=2)+
      geom_vline(xintercept=c(-0.02,0.02),linetype=2)+
      xlim(slope.L-0.1,slope.R+0.1)+
      geom_text_repel(force = 2,size=5)+theme_pubr()+ggtitle(Name)
  }else{
    p <- ggplot(datatoplot)+aes(slopes,-log10(ps),label=Label)+
      geom_point()+
      geom_hline(yintercept=-log10(pline),linetype=2)+
      geom_vline(xintercept=c(-0.02,0.02),linetype=2)+
      xlim(slope.L-0.1,slope.R+0.1)+
      geom_text_repel(force = 2,size=5)+theme_pubr()+ggtitle(Name)
  }
}
    

p1 <- PlotLinRegress_Vocano(LinOut.result.df)
p1


##


## Repeat regression with normalization ##

# Filter and normalize clone-level gene expression
LinOut.df.norm <- normalize_aggregated_expr(LinOut.df, scale = F)

# Aggregate gene expression 
boxplot(log(LinOut.df.norm[, c(14:44)]))

# Regression with normalized and scaled GE matrix
mode <- "poi" # "lm" "poi"
core <- as.numeric(8) # Number of cores, I.e. 8
LinOut.df_new <- LinOut.df.norm

## Run regression
if(mode=="lm"){
  print("run linear model")
  LinOut.result.norm <- Run_Lin_regression(LinOut.df_new,n.cores=core)
}else if(mode=="poi"){
  print("run poisson model")
  LinOut.result.norm <- Run_Lin_regression_poi(LinOut.df_new,n.cores=core)
}

LinOut.result.norm.df <- convert_results_to_df(LinOut.result.norm)

p2 <- PlotLinRegress_Vocano(LinOut.result.norm.df)
p2


##


## Regression based on condition
LinOut.result.cond <- Run_Lin_regression_poi.for_condition(LinOut.df_new)
LinOut.result.cond.df <- convert_results_to_df(LinOut.result.cond, x_cond = T)

p1 <- PlotLinRegress_Vocano(LinOut.result.cond.df$NT_NT)
print(p1)

p2 <- PlotLinRegress_Vocano(LinOut.result.cond.df$NT_AC)
print(p2)

p3 <- PlotLinRegress_Vocano(LinOut.result.cond.df$AC_AC)
print(p3)

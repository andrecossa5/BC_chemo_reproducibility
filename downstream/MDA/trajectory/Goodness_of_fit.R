
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


mode <- "poi"
core <- as.numeric(8) 
LinOut.df_rds <- fs::path(path_input_data, "agg_for_poisson.RDS") 
name <- "reg_out" 

## Read in the LinOut.gene.df file
print("Read in")
print(LinOut.df_rds)
LinOut.df<-readRDS(LinOut.df_rds)
LinOut.df.filt <- filter_genes(LinOut.df)

# hclust
df_to_clust <- LinOut.df.filt[, 14:dim(LinOut.df.filt)[2]]
#dist_matrix <- dist(t(df_to_clust))  
#hclust_result <- hclust(dist_matrix)  
#plot(hclust_result)

#num_clusters <- 4
#cluster_assignments <- cutree(hclust_result, k = num_clusters)

# avg. expression
clonsize.umi <- df_to_clust %>% rowSums()
df_to_clust.norm <-  1e6*df_to_clust / clonsize.umi
avg_expr <- apply(df_to_clust.norm, MARGIN = 2, FUN = mean)
plot(sort(avg_expr), type = "l")
hist(sort(avg_expr))

p25 <- quantile(avg_expr)[2]
p75 <- quantile(avg_expr)[4]


##


regress_factor <- "met_potential"
LinOut.df.filt$tot_UMIs <- LinOut.df.filt[, 14:dim(LinOut.df.filt)[2]] %>% rowSums()

# Fit model to 5 genes from set avg_expr < p25
low <-df_to_clust.norm[, avg_expr < p25]
low_to_check <- colnames(low)[1:5]

# Fit model to 5 genes from set avg_expr > p75
high <-df_to_clust.norm[, avg_expr > p75]
high_to_check <- colnames(high)[1:5]

# Fit model to 5 genes from set avg_expr > p25 & < p75
middle <-df_to_clust.norm[, avg_expr > p25 & avg_expr < p75]
mid_to_check <- colnames(middle)[1:5]

gene_sets <- list("low" = low_to_check, 
                  "mid" = mid_to_check, 
                  "high" = high_to_check)

# Verify goodness of models 
for(set in names(gene_sets)){
  gene_set <- gene_sets[[set]]
  
  cat("\n")
  print(paste0("--- Analyzing genes from set --- : ", set))
  for(gene in gene_set){
    cat("\n")
    print(gene)
    f <- as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)"))
    mdp <- glm(f, data=LinOut.df.filt,family=poisson(link="log"))
    mdnb <- MASS::glm.nb(f, data=LinOut.df.filt)
    
    ## Residuals
    # Poi
    residualsp <- residuals(mdp)
    p <- plot(mdp$fitted.values, residualsp, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted - Poi")
    print(p)
    # Nb
    residualsnb <- residuals(mdnb)
    p <- plot(mdnb$fitted.values, residualsnb, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted - Nb")
    print(p)
    
    ## Test 
    # Poi
    gof_p <- pchisq(sum(residualsp^2), df = df.residual(mdp), lower.tail = FALSE)
    print("Goodness of fit - Poi:")
    print(gof_p)
    # Nb
    gof_nb <- pchisq(sum(residualsnb^2), df = df.residual(mdnb), lower.tail = FALSE)
    print("Goodness of fit - Nb:")
    print(gof_nb)
  }
}


# The p-value obtained from pchisq represents the probability of observing a chi-squared value as extreme as or more extreme 
# than the calculated chi-squared statistic, given the null hypothesis that the model fits the data well.
# If the p-value is small (typically less than a chosen significance level, often 0.05), it suggests that the observed 
# chi-squared statistic is unlikely to have occurred by chance under the assumption of good model fit. In other words, 
# there's evidence to reject the null hypothesis, indicating poor model fit.



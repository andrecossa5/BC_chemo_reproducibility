
### Trajectory analysis on MDA clones ###

# Functions at: https://github.com/chenweng1991/redeemR/blob/master/R/BuidTree.R



#' Tomerge_v2
#'
#' This function is to quickly merge two dataframe by rownames, but can choose to leave A or B all information
#' @param A  dataframe A
#' @param B  dataframe B
#' @return  return a data frame with merged information
#' @export
#' @examples
#' Gettopgenes(XX.ob,number)
Tomerge_v2<-function(A,B,leavex=T,leavey=F){
  mergeAB<-merge(A,B,by="row.names",all.x=leavex,all.y=leavey)
  row.names(mergeAB)<-mergeAB[,1]
  mergeAB<-mergeAB[,-1]
  return(mergeAB)
}


#' MakeDF4Regress  
#' Define a function to make two dataframe for regression analysis
#' This function was developed based on HSC_multiome_Het_2.ipynb
#' @param multiome_wrapper This outject should includes all and more than HSCs cells in redeemR
#' @param redeemR  scredeemR object for HSC
#' @param progeny_np run via ProgenyMapping_np
#' @param assay SCT for expression, ATAC for ATAC
#' @param useNPimputation default is T, use all cells called by network propagation, inaddition to the top cells in redeemR
#' @param maxcloneUMI default is 10, Only include genes, in the max clone the expression greater than 10
#' @return list(mtx.clone=mtx.clone,mtx.clone.norm.scale=mtx.clone.norm.scale)
#' @export
MakeDF4Regress<-function(multiome_wrapper=Donor04_HSC_Multiome_wrapper,
                         redeemR=DN4_stemcell_redeemR.seed.sensitive,
                         progeny_np=DN4_HSC_LSI_progeny_np,
                         assay="SCT", useNPimputation=T, maxcloneUMI=10){
  
  ## Get mtx with counts (should be raw)
  mtx<-multiome_wrapper@assays[[assay]]@counts %>% t %>% as.matrix
  row.names(mtx)<-gsub("Cell","",row.names(mtx))
  
  ## Compute CloneInfo (Cell|npClone)
  if(useNPimputation){
    CloneInfo<-subset(progeny_np$ALLmeta.npClone,STD.CellType=="HSC") %>% .[,c("Cell","npClone")] %>% mutate(Cell_RNA=Translate_simple_ATAC2RNA(Cell)) 
  }else{
    CloneInfo<-redeemR@CellMeta[,c("Cell","Clone_merge")]%>% mutate(Cell_RNA=Translate_simple_ATAC2RNA(Cell)) %>% rename(npClone="Clone_merge")   
  }
  CloneInfo<-CloneInfo %>%tibble::remove_rownames() %>% tibble::column_to_rownames("Cell_RNA") %>% select(-Cell)
  
  ## Compute mtx.clone.norm(npClone|gene1|gene2|...)
  mtx.clone<-Tomerge_v2(mtx,CloneInfo) %>% group_by(npClone) %>% dplyr::summarise_all(sum) # aggregate counts by clone, by summing counts
  mtx.clone<-mtx.clone[!is.na(mtx.clone$npClone),]  ## Becasue some cells is NA when assigning clones
  clonesize.umi<-mtx.clone[,-1] %>% rowSums ## The order of the clonesize.umi is consistent with the rows of mtx.clone 
  mtx.clone.norm <- 1e6*mtx.clone[,-1]/clonesize.umi # Multiply for a factor and divide by 'clonesize'
  mtx.clone.norm.scale<-cbind(npClone=mtx.clone[,1],scale(mtx.clone.norm)) # Center by subtracting column means
  
  ## Only include genes with at least maxcloneUMI in max clone
  mtx.clone[,-1] %>% apply(.,2,max) -> maxcloneExpr
  features<-names(maxcloneExpr)[maxcloneExpr>=maxcloneUMI]
  mtx.clone<-mtx.clone[,c("npClone",features)] %>% mutate(npClone=as.character(npClone)) %>% cbind(TotalUMI=clonesize.umi,.)
  mtx.clone.norm.scale<-mtx.clone.norm.scale[,c("npClone",features)] %>% mutate(npClone=as.character(npClone))
  
  ## Merge in the lineage bias and output information
  Lin_Out<-progeny_np$datatoplot.scale %>% rename(npClone="Clone_merge")
  mtx.clone<-merge(Lin_Out,mtx.clone)
  mtx.clone.norm.scale<-merge(Lin_Out,mtx.clone.norm.scale)
  return(list(mtx.clone=mtx.clone,mtx.clone.norm.scale=mtx.clone.norm.scale))
}


#' Filter genes
#' 
#' This function is used to filter out genes with zero-expression in too many samples
#' 
#' @param 
#' 
#' @return 
filter_genes <- function(input_expr_df, perc_thr = 10){
  n_clones <- length(input_expr_df$GBC)
  n_clones_thr <- ceiling((n_clones / 100) * perc_thr)
  
  n_zeros <- colSums(input_expr_df[,14:dim(input_expr_df)[2]] == 0) 
  features <- names(n_zeros)[!n_zeros>n_clones_thr] # Retain only genes with non-zero expr in at least 90% of clones
  
  out_expr_df <- input_expr_df[, c(colnames(input_expr_df)[1:13], features)] 
  return(out_expr_df)
}


#' Convert to unique data.frame 
#' 
#' @param  
#' @return 
convert_results_to_df <- function(LinOut.result, x_cond = F){
  if(!x_cond){
    LinOut.result.df <- merge(LinOut.result$slopes, LinOut.result$ps, by = "row.names") %>%
      merge(., LinOut.result$qs, by.x = "Row.names", by.y = "row.names", all = T)
    colnames(LinOut.result.df) <- c("genes", "slopes", "ps", "qs")
    
    return(LinOut.result.df)    
  } else if(x_cond){
    for(exp_cond in names(LinOut.result)){
      LinOut.result.cond <- LinOut.result[[exp_cond]]
      LinOut.result.cond.df <- merge(LinOut.result.cond$slopes, LinOut.result.cond$ps, by = "row.names") %>%
        merge(., LinOut.result.cond$qs, by.x = "Row.names", by.y = "row.names", all = T)
      colnames(LinOut.result.cond.df) <- c("genes", "slopes", "ps", "qs")
      
      LinOut.result[[exp_cond]] <- LinOut.result.cond.df
    }
    return(LinOut.result)
  }
}


#' Pre-process aggregated gene-expression
#' 
#' This function is used to filter, normalize, and scale clone-level gene expression data 
#' 
#' @param  
#' @return 
normalize_aggregated_expr <- function(input_expr_df, verbose = F, log = T, scale = F){
  if(verbose){
    if(log){
      b <- boxplot(log(input_expr_df[, c(14:(14+30))]))
    } else {
      b <- boxplot(input_expr_df[, c(14:(14+30))])
    }
    print("Clone-level expression data before normalization")
    invisible(b)
  }
  
  # Filter out genes with zero-expression in too many clones
  input_expr_df.filt <- filter_genes(input_expr_df)
  
  # Normalize by total n. of counts x clone, round
  clonesize.umi <- input_expr_df.filt[,-(1:13)] %>% rowSums 
  input_expr_df.filt.norm <- 1e6*input_expr_df.filt[,-(1:13)]/clonesize.umi
  
  if(scale == T){
    # Scale, min-shift to have only positive values, and round to have integer counts (poisson)
    input_expr_df.filt.norm.scale <- scale(input_expr_df.filt.norm)
    input_expr_df.filt.norm.scale.min <- input_expr_df.filt.norm.scale + abs(min(input_expr_df.filt.norm.scale))
    input_expr_df.filt.norm.scale.min.round <- round(input_expr_df.filt.norm.scale.min)
    input_expr_df.final <- cbind(input_expr_df[,1:13], input_expr_df.filt.norm.scale.min.round)
  } else {
    # Only round CPMs 
    input_expr_df.filt.norm.round <- round(input_expr_df.filt.norm)
    input_expr_df.final <- cbind(input_expr_df[,1:13], input_expr_df.filt.norm.round)
  }
  
  if(verbose){
    if(log){
      b <- boxplot(log(input_expr_df.filt.norm.scale.min.round[, c(14:(14+30))]))
    } else {
      b <- boxplot(input_expr_df.filt.norm.scale.min.round[, c(14:(14+30))])
    }
    print("Clone-level expression data AFTER normalization")
    invisible(b)
  }
  
  return(input_expr_df.final)
}


#' Run_Lin_regression
#' 
#' Firstly used in HSC_multiome_Het_2.ipynb
#' @param LinOut  produced by MakeDF4Regress
#' @param regress_factor default is c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME")
#' @param n.cores  default is 8
#' @export
#' @import foreach doParallel
Run_Lin_regression<-function(LinOut,regress_factor=c("OutLevel.scale","OutLevel_NPadj.scale","Lym","Mye","MK","ME"),n.cores=8){
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  colnames(LinOut$mtx.clone.norm.scale)<-gsub("-","_",colnames(LinOut$mtx.clone.norm.scale))
  colnames(LinOut$mtx.clone.norm.scale)<-gsub("/","_",colnames(LinOut$mtx.clone.norm.scale))
  genes=colnames(LinOut$mtx.clone.norm.scale)[8:ncol(LinOut$mtx.clone.norm.scale)]
  my.cluster <- parallel::makeCluster(n.cores)
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){
      f<-as.formula(paste(gene,"~",regress_factor)[i])
      md<-lm(f,data=LinOut$mtx.clone.norm.scale)
      md.summary<-summary(md)
      slope<-md.summary$coefficients[2,1]
      p<-md.summary$coefficients[2,4]
      slopes<-c(slopes,slope)
      ps<-c(ps,p)
    }
    return(list(Gene=gene,ps=ps,slopes=slopes))
  }
  parallel::stopCluster(cl = my.cluster)
  genes<-sapply(res,function(x){x[[1]]})
  ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame
  slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame
  genes<-gsub("_","-",genes)
  row.names(ps)<-genes
  row.names(slopes)<-genes
  colnames(ps)<-regress_factor
  colnames(slopes)<-regress_factor
  qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame  
  return(list(ps=ps,qs=qs,slopes=slopes))
}



#' Run_Lin_regression_poi
#' 
#' Run Poisson regression on clone-level gene expression data
#' 
#' @param LinOut input data.frame with clone-level gene expression data plus metastatic potential
#' @param regress_factor  
#' @param n.cores  =8
#' @param qval compute q-values additionally to p-values. Sometimes q-values computation is not successful. Keep p-values then.
#' @param tot_UMIs include total n. of UMIs x clone instead of median. This allows to have a functioning poi model.
#' @export
#' @import foreach doParallel 
Run_Lin_regression_poi<-function(LinOut, 
                                 regress_factor = c("met_potential"),
                                 n.cores = 8, 
                                 qval = F, 
                                 tot_UMIs = F){ 
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  
  genes=colnames(LinOut)[14:ncol(LinOut)]
  
  # To allow the model to correct for clone-size, compute total n. of counts x clone
  if(tot_UMIs == T){
    clonesize.umi <- LinOut[,14:dim(LinOut)[2]] %>% rowSums 
    LinOut <- LinOut %>% mutate(tot_UMIs = clonesize.umi) %>% 
      relocate(., tot_UMIs, .after = nUMIs)
  }
  
  # used to set up a parallel computing environment
  my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
  
  res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
      if(tot_UMIs == T){
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
      }else{
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
      }
      md<-glm(f,data=LinOut,family=poisson(link="log"))
      md.summary<-summary(md)
      slope<-md.summary$coefficients[2,1] # Save slope from fitted model
      p<-md.summary$coefficients[2,4] # Save p-value
      slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
      ps<-c(ps,p) # Store each p-value
    }
    return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
  }
  parallel::stopCluster(cl = my.cluster)
  
  if(length(regress_factor) > 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
  } else if(length(regress_factor) == 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
  }
  
  row.names(ps)<-genes
  row.names(slopes)<-genes
  colnames(ps)<-regress_factor
  colnames(slopes)<-regress_factor
  
  # Save either q-values of p-values based on choice 
  if(qval == F){
    return(list(ps=ps,slopes=slopes))  
  }else{
    qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame
    return(list(ps=ps,qs=qs,slopes=slopes))
  }
}




#' Run_Lin_regression_poi
#' 
#' This function was developed based on 
#' @param LinOut produced by MakeDF4Regress
#' @param regress_factor  default is 
#' @param n.cores  =8
#' @export
#' @import foreach doParallel 
Run_Lin_regression_poi.for_condition<-function(LinOut, 
                                 regress_factor = c("met_potential"),
                                 n.cores = 8, 
                                 qval = F, 
                                 tot_UMIs = T){
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  
  
  # Split data.frame based on experimental condition
  LinOut <- LinOut %>% mutate(exp_condition = str_sub(LinOut$dataset, end = -3)) %>% 
    relocate(., exp_condition, .after = dataset) 
  LinOut_split <- split(LinOut, LinOut$exp_condition)
  
  res_x_cond <- list()
  for(exp_cond in names(LinOut_split)){
    # Run regression based on experimental condition
    LinOut_split_cond <- LinOut_split[[exp_cond]]

    genes=colnames(LinOut_split_cond)[15:ncol(LinOut_split_cond)]
    
    # To allow the model to correct for clone-size, compute total n. of counts x clone
    if(tot_UMIs == T){
      clonesize.umi <- LinOut_split_cond[,15:dim(LinOut_split_cond)[2]] %>% rowSums 
      LinOut_split_cond <- LinOut_split_cond %>% mutate(tot_UMIs = clonesize.umi) %>% 
        relocate(., tot_UMIs, .after = nUMIs)
    }
    
    # used to set up a parallel computing environment
    my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
    print(my.cluster)
    doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
    
    res<-foreach(gene=genes[1:length(genes)]) %dopar%{
      ps<-c()
      slopes<-c()
      for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
        if(tot_UMIs == T){
          f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
        }else{
          f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
        }
        md<-glm(f,data=LinOut_split_cond,family=poisson(link="log"))
        md.summary<-summary(md)
        slope<-md.summary$coefficients[2,1] # Save slope from fitted model
        p<-md.summary$coefficients[2,4] # Save p-value
        slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
        ps<-c(ps,p) # Store each p-value
      }
      return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
    }
    parallel::stopCluster(cl = my.cluster)
    
    if(length(regress_factor) > 1){
      genes<-sapply(res,function(x){x[[1]]})
      ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
      slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
    } else if(length(regress_factor) == 1){
      genes<-sapply(res,function(x){x[[1]]})
      ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
      slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
    }
    
    row.names(ps)<-genes
    row.names(slopes)<-genes
    colnames(ps)<-regress_factor
    colnames(slopes)<-regress_factor
    
    # Save either q-values of p-values based on choice 
    if(qval == F){
      res_x_cond[[exp_cond]] <- list(ps=ps,slopes=slopes)
    }else{
      qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame
      res_x_cond[[exp_cond]] <- list(ps=ps,qs=qs,slopes=slopes)
    }
  }
  return(res_x_cond)
}


#' Run_Lin_regression_poi
#' 
#' This function was developed based on 
#' @param LinOut 
#' @param regress_factor  default is 
#' @param n.cores  =8
#' @export
#' @import foreach doParallel 
Run_Lin_regression_nb.for_condition<-function(LinOut, 
                                               regress_factor = c("met_potential"),
                                               n.cores = 8, 
                                               qval = F, 
                                               tot_UMIs = T){
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  
  
  # Split data.frame based on experimental condition
  LinOut <- LinOut %>% mutate(exp_condition = str_sub(LinOut$dataset, end = -3)) %>% 
    relocate(., exp_condition, .after = dataset) 
  LinOut_split <- split(LinOut, LinOut$exp_condition)
  
  res_x_cond <- list()
  for(exp_cond in names(LinOut_split)){
    # Run regression based on experimental condition
    LinOut_split_cond <- LinOut_split[[exp_cond]]
    
    genes=colnames(LinOut_split_cond)[15:ncol(LinOut_split_cond)]
    
    # To allow the model to correct for clone-size, compute total n. of counts x clone
    if(tot_UMIs == T){
      clonesize.umi <- LinOut_split_cond[,15:dim(LinOut_split_cond)[2]] %>% rowSums 
      LinOut_split_cond <- LinOut_split_cond %>% mutate(tot_UMIs = clonesize.umi) %>% 
        relocate(., tot_UMIs, .after = nUMIs)
    }
    
    # used to set up a parallel computing environment
    my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
    print(my.cluster)
    doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
    
    res<-foreach(gene=genes[1:length(genes)]) %dopar%{
      ps<-c()
      slopes<-c()
      for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
        if(tot_UMIs == T){
          f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
        }else{
          f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
        }
        md<-MASS::glm.nb(f,data=LinOut_split_cond)
        md.summary<-summary(md)
        slope<-md.summary$coefficients[2,1] # Save slope from fitted model
        p<-md.summary$coefficients[2,4] # Save p-value
        slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
        ps<-c(ps,p) # Store each p-value
      }
      return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
    }
    parallel::stopCluster(cl = my.cluster)
    
    if(length(regress_factor) > 1){
      genes<-sapply(res,function(x){x[[1]]})
      ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
      slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
    } else if(length(regress_factor) == 1){
      genes<-sapply(res,function(x){x[[1]]})
      ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
      slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
    }
    
    row.names(ps)<-genes
    row.names(slopes)<-genes
    colnames(ps)<-regress_factor
    colnames(slopes)<-regress_factor
    
    # Save either q-values of p-values based on choice 
    if(qval == F){
      res_x_cond[[exp_cond]] <- list(ps=ps,slopes=slopes)
    }else{
      qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame
      res_x_cond[[exp_cond]] <- list(ps=ps,qs=qs,slopes=slopes)
    }
  }
  return(res_x_cond)
}



#' Run_Lin_regression_nb
#' 
#' Run Negative Binomial regression on clone-level gene expression data
#' Negative Binomial reg. is more suitable then Poisson reg. when overdispersion is present
#' 
#' @param LinOut input data.frame with clone-level gene expression data plus metastatic potential
#' @param regress_factor  
#' @param n.cores  =8
#' @param qval compute q-values additionally to p-values. Sometimes q-values computation is not successful. Keep p-values then.
#' @param tot_UMIs include total n. of UMIs x clone instead of median. This allows to have a functioning poi model.
#' @export
#' @import foreach doParallel 
Run_Lin_regression_nb<-function(LinOut, 
                                regress_factor = c("met_potential"),
                                n.cores = 8,
                                qval = F, 
                                tot_UMIs = T){
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(qvalue)
  library(MASS)
  
  genes=colnames(LinOut)[14:ncol(LinOut)]
  
  # To allow the model to correct for clone-size, compute total n. of counts x clone
  if(tot_UMIs == T){
    clonesize.umi <- LinOut[,14:dim(LinOut)[2]] %>% rowSums 
    LinOut <- LinOut %>% mutate(tot_UMIs = clonesize.umi) %>% 
      relocate(., tot_UMIs, .after = nUMIs)
  }
  
  # used to set up a parallel computing environment
  my.cluster <- parallel::makeCluster(n.cores) # This command creates a cluster of worker processes for parallel computation. n.cores specifies the number of CPU cores you want to use for parallel computation.
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster) # This command registers the parallel backend for the foreach package, allowing you to execute foreach loops in parallel.
  
  res<-foreach(gene=genes[1:length(genes)]) %dopar%{
    ps<-c()
    slopes<-c()
    for(i in 1:length(regress_factor)){ # For each factor for which we want to perform poi regression (i.e. for us, only clonal fitness)
      if(tot_UMIs == T){
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(tot_UMIs)")[i])
      }else{
        f<-as.formula(paste(gene,"~",regress_factor,"+ log(nUMIs)")[i])
      }
      md<-MASS::glm.nb(f,data=LinOut)
      md.summary<-summary(md)
      slope<-md.summary$coefficients[2,1] # Save slope from fitted model
      p<-md.summary$coefficients[2,4] # Save p-value
      slopes<-c(slopes,slope) # Store each slope, of model fitted on gene ~ each regress_factor
      ps<-c(ps,p) # Store each p-value
    }
    return(list(Gene=gene,ps=ps,slopes=slopes)) # res will be a list of lists (?)
  }
  parallel::stopCluster(cl = my.cluster)
  
  if(length(regress_factor) > 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% t %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% t %>% as.data.frame    
  } else if(length(regress_factor) == 1){
    genes<-sapply(res,function(x){x[[1]]})
    ps<-sapply(res,function(x){x[[2]]}) %>% as.data.frame 
    slopes<-sapply(res,function(x){x[[3]]}) %>% as.data.frame
  }
  
  row.names(ps)<-genes
  row.names(slopes)<-genes
  colnames(ps)<-regress_factor
  colnames(slopes)<-regress_factor
  
  # Save either q-values of p-values based on choice 
  if(qval == F){
    return(list(ps=ps,slopes=slopes))  
  }else{
    qs<-apply(ps,2,function(x){qvalue(x)$qvalues})  %>% as.data.frame
    return(list(ps=ps,qs=qs,slopes=slopes))
  }
}



#' Volcano plot of results
#' 
#' Plot a Volcano with q-value ~ slopes of genes resulting from poi regression
#' 
#' @param 
#' @return 
PlotLinRegress_Vocano <- function(LinOut.result.df, slot = "qs", pline = 0.001, qline = 0.05){
  library(ggpubr)
  library(purrr)
  
  datatoplot <- LinOut.result.df
  ifelse(slot == "qs", datatoplot <- datatoplot %>% mutate(score = -log(qs)*abs(slopes)), 
         datatoplot <- datatoplot %>% mutate(score = -log(ps)*abs(slopes)))
  
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
    # Cap p-values to min p-value to avoid having zeros (leading to -log10(0) = Inf)
    datatoplot$ps[datatoplot$ps == 0] <- min(datatoplot$ps[datatoplot$ps != 0])
    pval.top <- -log10(datatoplot$ps) %>% max()
    p <- ggplot(datatoplot)+aes(slopes,-log10(ps),label=Label)+
      geom_point()+
      geom_hline(yintercept=-log10(pline),linetype=2)+
      geom_vline(xintercept=c(-0.02,0.02),linetype=2)+
      xlim(slope.L-0.1,slope.R+0.1)+
      ylim(0, pval.top+0.1)+
      geom_text_repel(force = 2,size=5)+theme_pubr()+ggtitle(Name)
  }
}



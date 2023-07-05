
#' Calculate summary statistics for simulated and real bulk expression
#'
#' @param bulkSimulator_res a list of simulated bulk object from bulkSimulator() function
#' @param baseline_expr a real bulk expression matrix to compare against, with genes in rows and samples in columns.
#' @param var.genes genes used to calculate gene-wise correlation and sample-wise correlation. If NULL, will automatically
#'    choose top most variable genes from baseline_expr as var.genes
#' @param top_n_variable a numeric value indicating the number of top variable genes to be selected. Required when var.genes = NULL
#' @param n.core number of cores to use for parallel programming. Default = 1.
#'
#' @return a list of summary statistics for the simulated and real bulk expression. The list includes a summary statistics table that consist of
#'    the mean, variance, coefficient of variation (CV), and Median Absolute Deviation (MAD) for all genes.
#'    It also provides gene-wise and sample-wise pairwise correlations for var.genes.
#' @export
simulator_Statistics <-function(bulkSimulator_res,
                                baseline_expr = NULL,
                                var.genes = NULL,
                                top_n_variable = 300,
                                n.core = 1){
  if(max(baseline_expr)<100){
    stop('baseline expression should not be log-transformed')
  }

  baseline_expr = as.matrix(baseline_expr)

  row_names_list <- lapply(bulkSimulator_res, function(x) rownames(x[[1]]))
  common_row_names <- Reduce(intersect, row_names_list)

  cms = intersect(rownames(baseline_expr),common_row_names)
  if(length(cms)<top_n_variable){
    stop('Too few intersection genes found between the baseline expression and the simulated expression.
         Please check the gene symbol format used in the baseline expression or consider using a different baseline expression.')
  }

  baseline_expr = baseline_expr[cms,]

  bulkSimulator_methods = names(bulkSimulator_res)
  for(method in bulkSimulator_methods){
    bulkSimulator_res[[method]]$simulated_bulk = bulkSimulator_res[[method]]$simulated_bulk[cms,]
  }

  baseline_expr = apply(baseline_expr,2,function(x)((x/sum(x))*1e+06))

  bulk_list =  lapply(bulkSimulator_res, function(component) {
    component$simulated_bulk  # Extract the first element of each sub-component
  })

  bulk_list$baseline = baseline_expr

  message('Calculating summary statistics: mean, variance, cv and MAD')
  quick_cv<-function(v){
    sd(v)/mean(v)
  }

  quick_bulk_summary<-function(bulk_expr){
    df=data.frame(mean=apply(log2(bulk_expr+1),1,mean),
                  var=apply(log2(bulk_expr+1),1,var),
                  cv=apply(log2(bulk_expr+1),1,quick_cv),
                  MAD = apply(log2(bulk_expr+1),1,stats::mad),
                  frac_of_zero=apply(bulk_expr,1,function(x)sum(x==0))/ncol(bulk_expr),
                  gene=rownames(bulk_expr))
    return(df)
  }
  pboptions(type = "txt", style = 3, char = "=")
  summary_stats = do.call(rbind,pblapply(bulk_list,quick_bulk_summary,cl = n.core))
  summary_stats$class=sub("\\..*","",rownames(summary_stats))
  rownames(summary_stats)=NULL

  if(is.null(var.genes)){
    # select top most variable genes for pairwise correlation statistics
    summary_stats_baseline = summary_stats[summary_stats$class == 'baseline',]
    summary_stats_baseline = summary_stats_baseline[order(summary_stats_baseline$var,decreasing = T),]
    var.genes = summary_stats_baseline$gene[1:top_n_variable]
  }else{
    var.genes = var.genes[var.genes %in% cms]
  }

  message('Calculating gene correlations')
  quick_gene_cor = function(bulk_expr){
    return(cor(t(bulk_expr[var.genes,]),use = 'pairwise.complete.obs'))
  }
  pboptions(type = "txt", style = 3, char = "=")
  gene_correlation = pblapply(bulk_list,quick_gene_cor,cl = n.core)

  message('Calculating sample correlations')
  quick_sample_cor = function(bulk_expr){
    return(cor(bulk_expr[var.genes,],use = 'pairwise.complete.obs'))
  }
  pboptions(type = "txt", style = 3, char = "=")
  sample_correlation = pblapply(bulk_list,quick_sample_cor,cl = n.core)

  return(list(summary_stats = summary_stats,
              gene_correlation = gene_correlation,
              sample_correlation = sample_correlation))
}

#' Summarize all simulated fractions in simulated bulk object
#'
#' @param bulkSimulator_res a list of simulated bulk object from bulkSimulator() function
#' @param baseline_frac fraction matrix in real bulk expression with cell types in columns and samples in rows. Set to NULL if not available
#'
#' @return a dataframe that documents all the simulated fractions within the simulated bulk object. This dataframe is particularly useful for downstream visualization purposes
#' @export
simulator_fracStatistics = function(bulkSimulator_res, baseline_frac = NULL){
  frac_long = function(x){
    frac = x$simulated_frac %>% as.data.frame() %>% gather(cell_types,frac)
    return(frac)
  }

  f = do.call(rbind,lapply(bulkSimulator_res,frac_long))
  f$class=sub("\\..*","",rownames(f))
  rownames(f)=NULL

  if(!is.null(baseline_frac)){
    f_baseline = baseline_frac %>% as.data.frame() %>% gather(cell_types,frac) %>% mutate(class = 'baseline')
    f = rbind(f,f_baseline)
  }
  return(f)
}


#' Average gene-set CV calculation for simulated and real bulk expression
#'
#' @param simulator_Statistics_res a list of summary statistics retunred from simulator_Statistics() function
#' @param genelist a user-provided gene list. If set to NULL, will automatically retrieve gene sets collection from msigdbr
#' @param msigdbr_species msigdbr() species parameter. Required when genelist = NULL. Default = 'Homo sapiens'
#' @param msigdbr_category msigdbr() category parameter. Required when genelist = NULL. Default = 'H'
#' @param msigdbr_subcategory msigdbr() subcategory parameter.
#'
#' @return a dataframe of average CV for each gene set
#' @export
simulator_averageCV = function(simulator_Statistics_res, genelist = NULL, msigdbr_species = 'Homo sapiens',
                                msigdbr_category = 'H',msigdbr_subcategory = NULL){
  if(is.null(genelist)){
    message(paste('genelist not provided, will retrieve gene sets from msigdbr collections (default category: hallmark gene sets)'))
    if(!require(msigdbr,quietly = T)){
      stop('Suggested package "msigdbr" not installed')
    }else{
      require(msigdbr) %>% suppressMessages()
      retrieved_df = msigdbr(msigdbr_species, msigdbr_category, msigdbr_subcategory)
      genelist = msigdbr_to_list(retrieved_df)
    }
  }
  averageCV = genelist_cv(genelist,simulator_Statistics_res)
  return(averageCV)
}



#' Comparison of variance: simulated bulk vs real bulk
#'
#' @param simulator_Statistics_res a list of summary statistics returned from simulator_Statistics() function
#' @param genes_to_evalu genes to compare. If set to NULL, will compare all the genes in simulator_Statistics_res
#' @param variance_type a string indicating type of variance to compare. Available options including "cv" and "MAD". Deafult = "cv"
#' @param nrow_panel nrow of the ggplot panels
#' @param axis_limit axis_limit in the variance comparison plot
#' @param desired_order a character vector indicating desired order to organize simulated bulk in the figure.
#'    The default order is specified as c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC')
#' @param ... additional arguments to pass to ggpointdensity()
#'
#' @return a panel of ggpointdensity plots comparing variance in each simulated bulk with real bulk expression
#' @export
plot_variance_comparison = function(simulator_Statistics_res, genes_to_evalu = NULL, variance_type = 'cv',
                                    nrow_panel = 1,axis_limit = 5,
                                    desired_order = c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
                                    ...){

  summary_stats = simulator_Statistics_res$summary_stats
  if(!is.null(genes_to_evalu)){
    summary_stats = summary_stats[summary_stats$gene %in% genes_to_evalu,]
  }

  if(!require(ggpointdensity,quietly = T)){
    stop('Suggested package "ggpointdensity" not installed')
  }else{
    require(ggpointdensity)
  }

  summary_stats_simulatedExpr = summary_stats[summary_stats$class!= 'baseline',]
  summary_stats_baselineExpr = summary_stats[summary_stats$class == 'baseline',]

  df = summary_stats_simulatedExpr[,c('class',variance_type)]
  colnames(df)[2] = 'variance_simulatedExpr'
  df$variance_baselineExpr = rep(summary_stats_baselineExpr[,variance_type],length(unique(summary_stats_simulatedExpr$class)))

  df$class = factor(df$class,levels = desired_order)

  ggplot(df,aes(variance_baselineExpr,variance_simulatedExpr))+
    geom_pointdensity(...)+
    geom_abline(intercept = 0, slope = 1,linetype='dashed',color='red')+
    scale_x_continuous(limits = c(0, axis_limit), breaks = seq(0, axis_limit, 1)) +
    scale_y_continuous(limits = c(0, axis_limit), breaks = seq(0, axis_limit, 1))+
    facet_wrap(~class, nrow = nrow_panel)+
    theme_bw()+
    xlab(paste(variance_type,'in baseline expression'))+
    ylab(paste(variance_type,'in simulated bulk expression'))
}


#' Plot pairwise gene correlation in simulated bulk and real bulk expression
#'
#' @param simulator_Statistics_res a list of summary statistics returned from simulator_Statistics() function
#' @param nrow_panel nrow of the ggplot panels
#' @param desired_order a character vector indicating desired order to organize simulated bulk in the figure.
#'    The default order is specified as c('baseline','homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC')
#' @param ... additional arguments to pass to geom_tile()
#'
#' @return a panel of heatmaps showing pairwise gene correlations
#' @export
plot_PairwiseGeneCorrelation = function(simulator_Statistics_res, nrow_panel = 2,
                                        desired_order = c('baseline','homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
                                        ...){

  gene_cor = simulator_Statistics_res$gene_correlation
  distance_matrix <- 1 - gene_cor$baseline
  hclust_result <- hclust(dist(distance_matrix), method = "complete")
  hclust_order = hclust_result$order

  cor_mat_list <- lapply(gene_cor, function(x){
    x[hclust_result$order, hclust_result$order]})
  cor_mat_list <- lapply(cor_mat_list, reshape2::melt)

  cor_dat = do.call(rbind,lapply(cor_mat_list,function(x)x))
  cor_dat$class = sub("\\..*","",rownames(cor_dat))
  rownames(cor_dat) = NULL

  cor_dat$class = factor(cor_dat$class,levels = desired_order)

  ggplot(cor_dat, aes(Var2, Var1, fill = value))+
    facet_wrap(~class, nrow = nrow_panel) +
    geom_tile(...)+
    scale_fill_gradient2(low = "#477bb7", high = "#f73025", mid = "white",midpoint = 0, limit = c(-1,1), space = "Lab",name=" ") +

    theme(strip.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    xlab("") + ylab("") + coord_fixed()

}


#' Plot sample wise correlations in simulated bulk and real bulk
#'
#' @param simulator_Statistics_res a list of summary statistics returned from simulator_Statistics() function
#' @param desired_order a character vector indicating desired order to organize simulated bulk in the figure.
#'    The default order is specified as c('baseline','homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC')
#' @param ... additional arguments to pass to geom_violin()
#'
#' @return a panel of violin plot showing the distribution of sample-wise correlations
#' @export
plot_PairwiseSampleCorrelation = function(simulator_Statistics_res,
                                          desired_order = c('baseline','homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
                                          ...){
  sample_cor = simulator_Statistics_res$sample_correlation
  sample_cor = lapply(sample_cor,function(x){
    x = x[lower.tri(x,diag = F)]
  })

  cor_dat = data.frame(cor = unlist(sample_cor) %>% unname(),
                       class = rep(names(sample_cor),lengths(sample_cor)))

  cor_dat$class = factor(cor_dat$class,levels = desired_order)

  ggplot(cor_dat,aes(class,cor))+
    geom_violin(aes(color=class),...)+
    geom_boxplot(width=0.1,outlier.shape = NA,aes(color=class))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = 'none')+
    xlab('')+
    ylab('sample wise correlation')

}


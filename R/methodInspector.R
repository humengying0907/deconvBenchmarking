
#' Calculate summary statistics for simulated and real bulk expression
#'
#' @param bulkSimulator_res a list of simulated bulk object from bulkSimulator() function
#' @param baseline_expr a real bulk expression matrix to compare against, with genes in rows and samples in columns
#' @param var.genes genes used to calculate gene-wise correlation and sample-wise correlation. If NULL, will automatically
#'    choose top most variable genes from baseline_expr as var.genes
#' @param top_n_variable a numeric value indicating the number of top variable genes to be selected. Required when var.genes = NULL
#' @param n.core number of cores to use for parallel programming. Default = 1.
#'
#' @return a list of summary statistics for the simulated and real bulk expression. The list includes summary statistics such as
#'    mean, variance, coefficient of variation (CV), Median Absolute Deviation (MAD), as well as gene-wise and sample-wise pairwise correlations.
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
  cms = commonRows(baseline_expr,bulkSimulator_res[[1]]$simulated_bulk)
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


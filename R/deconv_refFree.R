
#' Reference free deconvolution with linseed
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param k number of cell types in bulk expression
#'
#' @return a matrix of deconvolution results with cell types in rows and samples in columns
#' @export
deconv_refFree_linseed = function(bulk_expr,k){
  require(linseed,quietly = T)
  l <- LinseedObject$new(bulk_expr)
  l$calculatePairwiseLinearity()
  l$calculateSpearmanCorrelation()
  l$calculateSignificanceLevel(100)
  l$filterDatasetByPval(0.01)

  l$setCellTypeNumber(k)
  l$project("full") # projecting full dataset
  l$project("filtered")
  l$smartSearchCorners(dataset="filtered", error="norm")

  l$deconvolveByEndpoints()
  return(l$proportions)
}


#' Reference free deconvolution with debCAM
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param k number of cell types in bulk expression
#' @param corner.strategy The method to find corners of convex hull. 1: minimum sum of margin-of-errors; 2: minimum sum of reconstruction errors. The default is 2.
#' @param dim.rdc Reduced data dimension; should be not less than maximum candidate K.
#' @param thres.low The lower bound of percentage of genes to keep for CAM with ranked norm. The value should be between 0 and 1. The default is 0.05.
#' @param thres.high 	The higher bound of percentage of genes to keep for CAM with ranked norm. The value should be between 0 and 1. The default is 0.95.
#' @param cluster.method The method to do clustering. The default "K-Means" will use kmeans. The alternative "apcluster" will use apclusterK-methods.
#' @param cluster.num The number of clusters; should be much larger than K. The default is 50.
#' @param MG.num.thres The clusters with the gene number smaller than MG.num.thres will be treated as outliers. The default is 20.
#' @param lof.thres Remove local outlier using lofactor. MG.num.thres is used as the number of neighbors in the calculation of the local outlier factors. The default value 0.02 will remove top 2% local outliers. Zero value will disable lof.
#' @param quickhull Perform quickhull to select clusters or not. The default is True.
#' @param quick.select The number of candidate corners kept after quickhull and SFFS greedy search. If Null, only quickhull is applied. The default is 20. If this value is larger than the number of candidate corners after quickhull, greedy search will also not be applied.
#' @param sample.weight Vector of sample weights. If NULL, all samples have the same weights. The length should be the same as sample numbers. All values should be positive.
#' @param appro3 	Estimate A and S matrix by approach 3 or not. Please see CAMASest for further information. The default is TRUE.
#' @param generalNMF 	If TRUE, the decomposed proportion matrix has no sum-to-one constraint for each row. The default is FALSE. TRUE value brings two changes: (1) Without assuming samples are normalized,
#'    the first principal component will not forced to be along c(1,1,..,1) but a standard PCA will be applied during preprocessing. (2) Without sum-to-one constraint for each row, the scale ambiguity of each column vector in proportion matrix will not be removed.
#' @param cores The number of system cores for parallel computing. If not provided, one core for each element in K will be invoked. Zero value will disable parallel computing.
#'
#' @return a matrix of deconvolution results with cell types in rows and samples in columns
#' @export
deconv_refFree_debCAM = function(bulk_expr, k, corner.strategy = 2, dim.rdc = 10,
                                 thres.low = 0.05, thres.high = 0.95, cluster.method = "K-Means", cluster.num = 50, MG.num.thres = 20,
                                 lof.thres = 0.02, quickhull = TRUE, quick.select = NULL,
                                 sample.weight = NULL, appro3 = TRUE, generalNMF = FALSE,
                                 cores = NULL){
  require(debCAM,quietly = T)

  rCAM=debCAM::CAM(bulk_expr, K = k, corner.strategy, dim.rdc,
                   thres.low, thres.high, cluster.method, cluster.num, MG.num.thres,
                   lof.thres, quickhull, quick.select,
                   sample.weight, appro3, generalNMF,
                   cores)
  rCAM_proportion=debCAM::Amat(rCAM, k) %>% t()
  colnames(rCAM_proportion)=colnames(bulk_expr)
  rownames(rCAM_proportion)=paste0('cellType_',seq(1,k))

  return(rCAM_proportion)
}


#' Reference free deconvolution
#'
#' @param methods a character vector indicating which methods to use. Use list_deconv_refFree() to check for available method names.
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param k number of cell types in bulk expression
#' @param corner.strategy The method to find corners of convex hull. 1: minimum sum of margin-of-errors; 2: minimum sum of reconstruction errors. The default is 2.
#' @param dim.rdc Reduced data dimension; should be not less than maximum candidate K.
#' @param thres.low The lower bound of percentage of genes to keep for CAM with ranked norm. The value should be between 0 and 1. The default is 0.05.
#' @param thres.high 	The higher bound of percentage of genes to keep for CAM with ranked norm. The value should be between 0 and 1. The default is 0.95.
#' @param cluster.method The method to do clustering. The default "K-Means" will use kmeans. The alternative "apcluster" will use apclusterK-methods.
#' @param cluster.num The number of clusters; should be much larger than K. The default is 50.
#' @param MG.num.thres The clusters with the gene number smaller than MG.num.thres will be treated as outliers. The default is 20.
#' @param lof.thres Remove local outlier using lofactor. MG.num.thres is used as the number of neighbors in the calculation of the local outlier factors. The default value 0.02 will remove top 2% local outliers. Zero value will disable lof.
#' @param quickhull Perform quickhull to select clusters or not. The default is True.
#' @param quick.select The number of candidate corners kept after quickhull and SFFS greedy search. If Null, only quickhull is applied. The default is 20. If this value is larger than the number of candidate corners after quickhull, greedy search will also not be applied.
#' @param sample.weight Vector of sample weights. If NULL, all samples have the same weights. The length should be the same as sample numbers. All values should be positive.
#' @param appro3 	Estimate A and S matrix by approach 3 or not. Please see CAMASest for further information. The default is TRUE.
#' @param generalNMF 	If TRUE, the decomposed proportion matrix has no sum-to-one constraint for each row. The default is FALSE. TRUE value brings two changes: (1) Without assuming samples are normalized,
#'    the first principal component will not forced to be along c(1,1,..,1) but a standard PCA will be applied during preprocessing. (2) Without sum-to-one constraint for each row, the scale ambiguity of each column vector in proportion matrix will not be removed.
#' @param cores The number of system cores for parallel computing. If not provided, one core for each element in K will be invoked. Zero value will disable parallel computing.
#'
#' @return a list of deconvolution results obtained using different reference free methods. The results are named using the format 'RefFree_' + method
#' @export
#'
#' @examples
#' \dontrun{
#' deconv_refFree(methods = c('linseed','debCAM'),
#'                bulk_expr = bulk_expr, k = 8)
#' }
deconv_refFree = function(methods,
                          bulk_expr, k,
                          corner.strategy = 2, dim.rdc = 10,
                          thres.low = 0.05, thres.high = 0.95, cluster.method = "K-Means", cluster.num = 50, MG.num.thres = 20,
                          lof.thres = 0.02, quickhull = TRUE, quick.select = NULL,
                          sample.weight = NULL, appro3 = TRUE, generalNMF = FALSE,
                          cores = NULL){

  l = list_deconv_refFree(show_description = T)
  refRefDeconv_list = list()

  for(method in methods){
    description = l$description[l$method == method]
    message(description)
    switch(method,
           linseed = {
             result = list(RefFree_linseed = deconv_refFree_linseed(bulk_expr,k))
           },
           debCAM = {
             result = tryCatch({
               list(RefFree_debCAM = deconv_refFree_debCAM(bulk_expr, k, corner.strategy, dim.rdc,
                                                           thres.low, thres.high, cluster.method, cluster.num, MG.num.thres,
                                                           lof.thres, quickhull, quick.select,
                                                           sample.weight, appro3, generalNMF,
                                                           cores))
             }, error = function(e) {
               message("Error in debCAM method: ", e$message)
               NULL
             })
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_deconv_regression() to check for available methods"))
             result <- NULL
           })
    refRefDeconv_list = c(refRefDeconv_list,result)

  }
  return(refRefDeconv_list)
}





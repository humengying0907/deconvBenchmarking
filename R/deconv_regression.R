
#' Regression-based deconvolution with nnls
#' @description Perform bulk deconvolution using the Non-Negative Least Squares (NNLS) algorithm with input signature matrices
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param sigMatrix_list a list of signature matrix
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input signature matrices. The results are named using the format
#'    'RefBased_nnls_' followed by the name of the respective input signature matrix.
#' @export
deconv_regression_nnls = function(bulk_expr,sigMatrix_list, n.core = 1){

  quick_nnls = function(ref){
    cms = commonRows(bulk_expr,ref)
    bulk_expr = bulk_expr[cms,]
    ref = ref[cms,]
    nnls_res =do.call(cbind,lapply(apply(bulk_expr,2,function(x) nnls::nnls(ref,x)), function(y) y$x))
    nnls_res = sweep(nnls_res,2,colSums(nnls_res),'/')
    rownames(nnls_res) <- colnames(ref)
    return(list(nnls_res))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(sigMatrix_list,quick_nnls, cl = n.core))
  names(out) = paste0('RefBased_nnls_',names(out))
  return(out)
}


#' Regression-based deconvolution with CIBERSORT
#' @description Perform bulk deconvolution using the cibersort algorithm with input signature matrices
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param sigMatrix_list a list of signature matrix
#' @param cibersort_path path to CIBERSORT.R
#' @param QN a logical variable determining whether to quantile normalize the input matrix. Default = F
#' @param skip_raw a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources.
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input signature matrices. The results are named using the format
#'    'RefBased_cibersort_' followed by the name of the respective input signature matrix.
#' @export
deconv_regression_cibersort = function(bulk_expr,sigMatrix_list, cibersort_path = NULL, QN = F, skip_raw = T, n.core = 1){
  if(is.null(cibersort_path)){
    stop('please provide path to cibersort R script to run cibersort')
  }

  if(skip_raw == T){
    sigMatrix_list = sigMatrix_list[names(sigMatrix_list)!='raw']
  }

  require(e1071,quietly = T)
  require(parallel,quietly = T)
  require(preprocessCore,quietly = T)
  source(paste0(cibersort_path,'/CIBERSORT.R'))

  quick_cbs = function(ref){
    cms = commonRows(bulk_expr,ref)
    bulk_expr = bulk_expr[cms,]
    ref = ref[cms,]

    cbs_res = CIBERSORT(sig_matrix = ref, mixture_file = bulk_expr,QN = QN)
    cbs_res = t(cbs_res[,1:(ncol(cbs_res)-3)])
    return(list(cbs_res))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(sigMatrix_list,quick_cbs, cl = n.core))
  names(out) = paste0('RefBased_cibersort_',names(out))
  return(out)
}


#' Regression-based deconvolution with MuSiC
#' @description Perform bulk deconvolution using MuSiC with genes selected from the input signature matrices as MuSiC input markers
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param sigMatrix_list a list of signature matrix
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param normalize MuSiC::music_prop() normalize argument. Default = F
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different marker genes from input signature matrices. The results are named using the format
#'    'RefBased_MuSiC_' followed by the name of the respective input signature matrix.
#' @export
deconv_regression_MuSiC = function(bulk_expr, sigMatrix_list, scExpr = NULL, scMeta = NULL,
                                   colnames_of_cellType = NA, colnames_of_sample = NA, normalize = F, n.core = 1){
  require(MuSiC,quietly = T)
  require(SingleCellExperiment,quietly = T)

  if(is.null(scExpr) | is.null(scMeta) | is.na(colnames_of_cellType) | is.na(colnames_of_sample)){
    stop('please provide required arguments to run MuSiC')
  }

  sc.sce = SingleCellExperiment::SingleCellExperiment(list(counts=scExpr),colData = scMeta)

  quick_MuSiC = function(ref){
    sig_genes = commonRows(bulk_expr,ref)
    res.list = music_prop(bulk.mtx = bulk_expr, sc.sce = sc.sce, clusters = colnames_of_cellType,
                           markers = sig_genes, normalize = normalize, samples = colnames_of_sample,
                           verbose = F)
    music_res = t(res.list$Est.prop.weighted)
    return(list(music_res))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(sigMatrix_list,quick_MuSiC, cl = n.core))
  names(out) = paste0('RefBased_MuSiC_',names(out))
  return(out)
}


#' Regression-based deconvolution with weighted robust linear regression
#' @description Perform bulk deconvolution using weighted robust linear regression (w-RLM) algorithm from LinDeconSeq package, with the input signature matrices
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param sigMatrix_list a list of signature matrix
#' @param skip_raw a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources.
#' @param weight LinDeconSeq::deconSeq() argument, a logical variable determining whether to weight gene or not. Default = T
#' @param intercept LinDeconSeq::deconSeq() argument, a logical variable determining whether to add intercept when using robust linear regression. Default = T
#' @param scale LinDeconSeq::deconSeq() argument, a logical variable determining whether to scale bulk gene expression or not. Default = T
#' @param QN LinDeconSeq::deconSeq() argument, a logical variable determining whether to quantile normalize bulk expression profile. Default = T
#' @param verbose LinDeconSeq::deconSeq() argument, a logical variable determining whether to print detailed information. Default = F
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input signature matrices. The results are named using the format
#'    'RefBased_wRLM_' followed by the name of the respective input signature matrix.
#' @export
deconv_regression_wRLM = function(bulk_expr,sigMatrix_list,
                                  skip_raw = TRUE,
                                  weight = TRUE,
                                  intercept = TRUE,
                                  scale = FALSE,
                                  QN = FALSE,
                                  verbose = FALSE,
                                  n.core = 1){
  require(LinDeconSeq,quietly = T)

  if(skip_raw == T){
    sigMatrix_list = sigMatrix_list[names(sigMatrix_list)!='raw']
  }

  quick_wRLM = function(ref){
    cms = commonRows(bulk_expr,ref)
    bulk_expr = bulk_expr[cms,]
    ref = ref[cms,]
    LinDeconSeq_res = LinDeconSeq::deconSeq(bulk_expr %>% as.data.frame(),
                                            ref %>% as.data.frame(),
                                            weight,intercept,scale,QN,verbose) %>% t()
    return(list(LinDeconSeq_res))
  }
  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(sigMatrix_list,quick_wRLM, cl = n.core))
  names(out) = paste0('RefBased_wRLM_',names(out))
  return(out)
}

#' Regression-based deconvolution with robust partial correlations
#' @description Perform bulk deconvolution using robust partial correlations (RPC) algorithm from EpiDISH package, with the input signature matrices
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param sigMatrix_list a list of signature matrix
#' @param skip_raw a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources.
#' @param maxit EpiDISH::epidish() argument, an integer indicating the limit of the number of IWLS iterations. Default = 100
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input signature matrices. The results are named using the format
#'    'RefBased_RPC_' followed by the name of the respective input signature matrix.
#' @export
deconv_regression_RPC = function(bulk_expr,sigMatrix_list, skip_raw = T, maxit = 100, n.core = 1){
  require(EpiDISH,quietly = T)

  if(skip_raw == T){
    sigMatrix_list = sigMatrix_list[names(sigMatrix_list)!='raw']
  }

  quick_RPC = function(ref){
    cms = commonRows(bulk_expr,ref)
    bulk_expr = bulk_expr[cms,]
    ref = ref[cms,]
    res.list = epidish(beta.m = bulk_expr, ref, method = "RPC", maxit)
    RPC_res = t(res.list$estF)
    return(list(RPC_res))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(sigMatrix_list,quick_RPC, cl = n.core))
  names(out) = paste0('RefBased_RPC_',names(out))
  return(out)
}


#' Regression-based deconvolution
#' @description Generate a list of regression-based deconvolution results with user-selected regression methods
#'
#' @param methods a character vector indicating which methods to use. Use list_deconv_regression() to check for available method names.
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param sigMatrix_list a list of signature matrix
#' @param cibersort_path path to CIBERSORT.R. This argument is required for 'cibersort' method
#' @param QN_cibersort a logical variable determining whether to quantile normalize the input matrix. Default = F. This argument is required for 'cibersort' method
#' @param skip_raw_cibersort a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources. This argument is required for 'cibersort' method
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns. This argument is required for 'MuSiC' method
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr. This argument is required for 'MuSiC' method
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta. This argument is required for 'MuSiC' method
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. This argument is required for 'MuSiC' method
#' @param normalize_MuSiC MuSiC::music_prop() normalize argument. Default = F. This argument is required for 'MuSiC' method
#' @param skip_raw_wRLM a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources. This argument is required for 'wRLM' method
#' @param weight_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to weight gene or not. Default = T. This argument is required for 'wRLM' method
#' @param intercept_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to add intercept when using robust linear regression. Default = T. This argument is required for 'wRLM' method
#' @param scale_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to scale bulk gene expression or not. Default = T. This argument is required for 'wRLM' method
#' @param QN_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to quantile normalize bulk expression profile. Default = T. This argument is required for 'wRLM' method
#' @param skip_raw_RPC a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources. This argument is required for 'RPC' method
#' @param maxit_RPC EpiDISH::epidish() argument, an integer indicating the limit of the number of IWLS iterations. Default = 100. This argument is required for RPC method
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input signature matrices and different regression methods. The results are named using the format
#'    'RefBased_' + regression method, followed by the name of the respective input signature matrix
#' @export
#'
#' @examples
#' \dontrun{
#' # create sigMatrix_list using marker genes derived from DE analysis limma and scran
#' sigMatrix_list = refMatrix(methods = c('limma','scran'), scExpr = scExpr, cell_type_labels = scMeta$cell_type)
#'
#' # explore other methods outside of the provided DE-based methods
#' help("pre_refMatrix_autogeneS")
#' help("pre_refMatrix_cibersortx")
#'
#' # include signature matrices generated from 'raw' methods for MuSiC, since MuSiC recommend to input scRNA profile without gene filtering
#' sigMatrix_list = refMatrix(methods = c('raw','limma','scran'),
#'                            scExpr = scExpr, cell_type_labels = scMeta$cell_type)
#'
#' deconv_regression(methods = c('nnls','cibersort','MuSiC','wRLM','RPC'),
#'                   bulk_expr = bulk_expr,sigMatrix_list = sigMatrix_list,
#'
#'                   # cibersort arguments
#'                   cibersort_path = 'scripts/',
#'
#'                   # MuSiC arguments
#'                   scExpr = scExpr,scMeta = scMeta,colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
#'
#'                   n.core = 2)
#' }
deconv_regression = function(methods,
                             bulk_expr,sigMatrix_list,
                             cibersort_path = NULL,
                             QN_cibersort = F, skip_raw_cibersort = T,
                             scExpr = NULL, scMeta = NULL, colnames_of_cellType = NA, colnames_of_sample = NA, normalize_MuSiC = F,
                             skip_raw_wRLM = TRUE, weight_wRLM = TRUE, intercept_wRLM = TRUE, scale_wRLM = FALSE, QN_wRLM = FALSE,
                             skip_raw_RPC = TRUE, maxit_RPC = 100,
                             n.core = 1){

  if('cibersort' %in% methods){
    if(is.null(cibersort_path)){
      stop('please provide path to cibersort R script to run cibersort')
    }
  }

  if('MuSiC' %in% methods){
    if(is.null(scExpr) | is.null(scMeta) | is.na(colnames_of_cellType) | is.na(colnames_of_sample)){
      stop('please provide required arguments to run MuSiC')
    }
  }

  l = list_deconv_regression(show_description = T)
  regDeconv_list = list()


  for (method in methods){
    description = l$description[l$method == method]
    message(description)

    switch(method,
           nnls = {
             result = deconv_regression_nnls(bulk_expr,sigMatrix_list, n.core)
           },
           cibersort = {
             result = deconv_regression_cibersort(bulk_expr,sigMatrix_list, cibersort_path, QN_cibersort, skip_raw_cibersort, n.core)
           },
           MuSiC = {
             result = deconv_regression_MuSiC(bulk_expr, sigMatrix_list, scExpr, scMeta, colnames_of_cellType, colnames_of_sample, normalize_MuSiC, n.core)
           },
           wRLM = {
             result = deconv_regression_wRLM(bulk_expr,sigMatrix_list, skip_raw_wRLM,weight_wRLM, intercept_wRLM, scale_wRLM, QN_wRLM, verbose = FALSE, n.core = n.core)
           },
           RPC ={
             result = deconv_regression_RPC(bulk_expr,sigMatrix_list, skip_raw_RPC, maxit_RPC, n.core)
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_deconv_regression() to check for available methods"))
             result <- NULL
           })
    regDeconv_list <- c(regDeconv_list, result)
  }

  return(regDeconv_list)
}

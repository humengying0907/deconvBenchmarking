
#' Marker-based deconvolution with first principal score
#' @description Perform bulk deconvolution by calculating the first principal score of each cell-type specific markers. Please note that the resulting score
#'    represents an abundance estimate and should not be interpreted as fractions.
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param marker_list a list of cell-type specific markers. Each list element is named by the corresponding method used to generate the markers.
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input markers. The results are named using the format
#'    'MarkerBased_firstPC_' followed by the name of the method used to generate the markers, such as 'limma'
#' @export
deconv_marker_firstPC = function(bulk_expr, marker_list, n.core = 1){

  quick_pca = function(markers_subList){
    firstPC = function(genes,bulk_expr){
      mat=bulk_expr[rownames(bulk_expr) %in% genes,]
      exp.PCA<-FactoMineR::PCA(t(as.matrix(mat)),graph = F,scale.unit = T)
      exp.PCA$ind$coord[,1]
    }
    firstPC_res = do.call(rbind,lapply(markers_subList, firstPC, bulk_expr))
    return(list(firstPC_res))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(marker_list,quick_pca, cl = n.core))
  names(out) = paste0('MarkerBased_firstPC_',names(out))
  return(out)
}

#' Marker-based deconvolution with GSVA
#' @description Calculate relative enrichment of cell-type specific markers using ssGSVA. Please note that the resulting score
#'    represents an abundance estimate and should not be interpreted as fractions.
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param marker_list a list of cell-type specific markers. Each list element is named by the corresponding method used to generate the markers.
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input markers. The results are named using the format
#'    'MarkerBased_gsva_' followed by the name of the method used to generate the markers, such as 'limma'
#' @export
deconv_markers_gsva = function(bulk_expr, marker_list, n.core = 1){
  require(GSVA)

  quick_gsva = function(markers_subList){
    gsva_scores=GSVA::gsva(bulk_expr,markers_subList,method='ssgsea',ssgsea.norm=F,verbose=F)
    return(list(gsva_scores))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(marker_list,quick_gsva, cl = n.core))
  names(out) = paste0('MarkerBased_gsva_',names(out))
  return(out)
}


#' Marker-based deconvolution with debCAM
#' @description Perform bulk deconvolution using AfromMarkers() function from debCAM package. The resulting scores have sum-to-one constraint and
#'    can be interpreted as fractions.
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param marker_list a list of cell-type specific markers. Each list element is named by the corresponding method used to generate the markers.
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input markers. The results are named using the format
#'    'MarkerBased_debCAM_' followed by the name of the method used to generate the markers, such as 'limma'
#' @export
deconv_marker_debCAM = function(bulk_expr, marker_list, n.core = 1){
  require(debCAM)
  quick_CAMmarker = function(markers_subList){
    CAMmarker_res = debCAM::AfromMarkers(bulk_expr,markers_subList)
    colnames(CAMmarker_res) = names(markers_subList)
    return(list(t(CAMmarker_res)))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(marker_list,quick_CAMmarker, cl = n.core))
  names(out) = paste0('MarkerBased_debCAM_',names(out))
  return(out)
}


#' Marker-based deconvolution with TOAST::MDeconv()
#' @description Perform bulk deconvolution using MDeconv() function from TOAST package. The resulting scores have sum-to-one constraint and
#'    can be interpreted as fractions.
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param marker_list a list of cell-type specific markers. Each list element is named by the corresponding method used to generate the markers.
#' @param alpha TOAST::MDeconv() argument. A vector including the prior mean for all cell types. Set to NULL to for TOAST/-P (partial reference-free deconvolution without prior)
#' @param sigma TOAST::MDeconv() argument. A vector including the prior standard deviation for all cell types.  Set to NULL to for TOAST/-P (partial reference-free deconvolution without prior)
#' @param epsilon TOAST::MDeconv() argument. A numeric variable to control the level of convergence. With a large epsilon, it converges quickly and the algorithm may not converge well.
#'    With a small epsilon, it converges slower and the algorithm may converge "too much". The default value is 1e-3, which we find is a reasonable threshold.
#' @param maxIter TOAST::MDeconv() argument. Number of maximum iterations.
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input markers. The results are named using the format
#'    'MarkerBased_TOAST_' followed by the name of the method used to generate the markers, such as 'limma'
#' @export
deconv_markers_TOAST = function(bulk_expr, marker_list, alpha = NULL, sigma = NULL,
                                epsilon = 0.001, maxIter = 1000,n.core = 1){
  require(TOAST)

  quick_PRF = function(markers_subList){
    PRF_res = TOAST::MDeconv(bulk_expr,markers_subList,alpha, sigma,epsilon,maxIter,verbose= FALSE)
    return(list(PRF_res$H))
  }

  pboptions(type = "txt", style = 3, char = "=")
  out = do.call(c,pblapply(marker_list,quick_PRF, cl = n.core))
  names(out) = paste0('MarkerBased_TOAST_',names(out))
  return(out)
}



#' Marker-based deconvolution
#' @description Generate a list of marker-based deconvolution results with user-selected methods
#'
#' @param methods a character vector indicating which methods to use. Use list_deconv_marker() to check for available method names.
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param marker_list a list of cell-type specific markers. Each list element is named by the corresponding method used to generate the markers.
#' @param alpha TOAST::MDeconv() argument. A vector including the prior mean for all cell types. Set to NULL to for TOAST/-P (partial reference-free deconvolution without prior)
#' @param sigma TOAST::MDeconv() argument. A vector including the prior standard deviation for all cell types.  Set to NULL to for TOAST/-P (partial reference-free deconvolution without prior)
#' @param epsilon TOAST::MDeconv() argument. A numeric variable to control the level of convergence. With a large epsilon, it converges quickly and the algorithm may not converge well.
#'    With a small epsilon, it converges slower and the algorithm may converge "too much". The default value is 1e-3, which we find is a reasonable threshold.
#' @param maxIter TOAST::MDeconv() argument. Number of maximum iterations.
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results obtained using different input markers and different marker-based methods. The results are named using the format
#'    'MarkerBased_' + marker-based method, followed by the name of the respective input markers
#' @export
#'
#' @examples
#' \dontrun{
#' # create marker list with refMarkers() function
#' marker_list = refMarkers(methods = c('limma','scran'),
#'               scExpr = scExpr,
#'               cell_type_labels = scMeta$cell_type)
#'
#' deconv_marker(methods = c('firstPC','gsva','debCAM','TOAST'),
#'               bulk_expr = bulk_expr,
#'               marker_list = marker_list)
#' }
deconv_marker = function(methods,
                         bulk_expr, marker_list,
                         alpha = NULL, sigma = NULL,
                         epsilon = 0.001, maxIter = 1000,n.core = 1){

  l = list_deconv_marker(show_description = T)

  markerDeconv_list = list()
  for (method in methods){
    description = l$description[l$method == method]
    message(description)

    switch(method,
           firstPC = {
             result = deconv_marker_firstPC(bulk_expr, marker_list, n.core)
           },
           gsva = {
             result = deconv_markers_gsva(bulk_expr, marker_list, n.core)
           },
           debCAM = {
             result = deconv_marker_debCAM(bulk_expr, marker_list, n.core)
           },
           TOAST = {
             result = deconv_markers_TOAST(bulk_expr, marker_list, alpha, sigma,
                                           epsilon, maxIter, n.core)
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_deconv_marker() to check for available methods"))
             result <- NULL
           })
    markerDeconv_list = c(markerDeconv_list,result)
  }
  return(markerDeconv_list)
}






#' List available bulk simulation strategies
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#' @return a dataframe of available bulk simulation strategies along with the corresponding function names to execute these methods.
#' @export
#'
#' @examples
#' list_bulkSimulator()
list_bulkSimulator = function(show_description = F){
  if(show_description==T){
    l = data.frame(method =  c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
                   function_to_call = c('bulkSimulator_homo','bulkSimulator_semi','bulkSimulator_heter','bulkSimulator_heter_sampleIDfree',
                                        'bulkSimulator_favilaco','bulkSimulator_immunedeconv','bulkSimulator_SCDC'),
                   suggested_packages = c('','','','scran','','immunedeconv','SCDC'),
                   description = c('Homogeneous bulk simulation',
                                   'Semi-heterogeneous bulk simulation',
                                   'Heterogeneous bulk simulation',
                                   'sampleID-independent heterogeneous bulk simulation',
                                   'Bulk simulation using Generator() function from Favilaco et al.',
                                   'Bulk simulation using make_bulk_eset() function from immunedeconv package',
                                   'Bulk simulation using generateBulk_norep() function from SCDC package'))
  }else if(show_description==F){
    l = data.frame(method =  c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
                   function_to_call = c('bulkSimulator_homo','bulkSimulator_semi','bulkSimulator_heter','bulkSimulator_heter_sampleIDfree',
                                        'bulkSimulator_favilaco','bulkSimulator_immunedeconv','bulkSimulator_SCDC'),
                   suggested_packages = c('','','','scran','','immunedeconv','SCDC'))

  }
  return(l)
}

#' List available fraction simulation strategies
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#'
#' @return a dataframe of available fraction simulation strategies along with the corresponding function names to execute these methods.
#' @export
#'
#' @examples
#' list_fracSimulator()
list_fracSimulator = function(show_description = F){
  if(show_description==T){
    l = data.frame(method = c('Beta','Dirichlet','favilaco','SCDC'),
                   function_to_call = c('fracSimulator_Beta','fracSimulator_Dirichlet',
                                        'fracSimulator_favilaco','fracSimulator_SCDC'),
                   suggested_packages = c('','','','SCDC'),
                   description = c('Simulate realistic cell-type fractions from beta distribution',
                                   'Simulate cell-type fractions from dirichlet distribution',
                                   'Simulate cell-type fractions using Generator() function from Favilaco et al.',
                                   'Simulate cell-type fractions using generateBulk_norep() function from SCDC package'))
  }else if(show_description==F){
    l = data.frame(method = c('Beta','Dirichlet','favilaco','SCDC'),
                   function_to_call = c('fracSimulator_Beta','fracSimulator_Dirichlet',
                                        'fracSimulator_favilaco','fracSimulator_SCDC'),
                   suggested_packages = c('','','','SCDC'))
  }
  return(l)
}


#' List available fraction to get cell-type markers from scRNA reference
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#'
#' @return a dataframe of available methods to get markers and the corresponding function names to execute these methods
#' @export
#'
#' @examples
#' list_refMarkers()
list_refMarkers = function(show_description = F){
  if(show_description==T){
    l = data.frame(method = c('limma','scran','sigMatrixList'),
                   function_to_call = c('refMarkers_limma','refMarkers_scran','refMarkers_sigMatrixList'),
                   suggested_packages = c('limma','BayesPrism',''),
                   description = c('Obtain cell-type specific markers with limma DE analysis',
                                   'Obtain cell-type specific markers with scran DE analysis',
                                   'Obtain cell-type specific markers from a list of signature matrices'))
  }else if(show_description ==F){
    l = data.frame(method = c('limma','scran','sigMatrixList'),
                   function_to_call = c('refMarkers_limma','refMarkers_scran','refMarkers_sigMatrixList'),
                   suggested_packages = c('limma','BayesPrism',''))
  }
  return(l)
}


#' List available methods to construct signature matrix from scRNA reference
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#'
#' @return a dataframe of available methods to construct signature matrix and the corresponding function names to execute these methods
#' @export
#'
#' @examples
#' list_refMarix()
list_refMarix = function(show_description = F){
  if(show_description ==T){
    l = data.frame(method = c('raw','limma','scran','markerList'),
                   function_to_call = c('refMatrix_raw','refMatirx_limma', 'refMatirx_scran', 'refMatirx_markerList'),
                   suggested_packages = c('','limma','BayesPrism',''),
                   description = c('Construct signature matrix using average expression across cell types',
                                   'Construct signature matrix using markers derived from limma DE analysis',
                                   'Construct signature matrix using markers derived from scran DE analysis',
                                   'Construct signature matrix with a pre-calculated marker list'))
  }else if(show_description ==F){
    l = data.frame(method = c('raw','limma','scran','markerList'),
                   function_to_call = c('refMatrix_raw','refMatirx_limma', 'refMatirx_scran', 'refMatirx_markerList'),
                   suggested_packages = c('','limma','BayesPrism',''))
  }
  return(l)
}

#' List available methods to perform regression-based deconvolution
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#'
#' @return a dataframe of available regression-based deconvolution methods along with the corresponding function names to execute these methods.
#' @export
#'
#' @examples
#' list_deconv_regression()
list_deconv_regression = function(show_description = F){
  if(show_description ==T){
    l = data.frame(method = c('nnls','cibersort','MuSiC','wRLM','RPC'),
                   function_to_call = c('deconv_regression_nnls',
                                        'deconv_regression_cibersort',
                                        'deconv_regression_MuSiC',
                                        'deconv_regression_wRLM',
                                        'deconv_RPC'),
                   suggested_packages = c('nnls','e1071/parallel/preprocessCore','MuSiC/SingleCellExperiment','LinDeconSeq','EpiDISH'),
                   description = c('Regression-based deconvolution with Non-Negative Least Squares',
                                   'Regression-based deconvolution with CIBERSORT',
                                   'Regression-based deconvolution with MuSiC',
                                   'Regression-based deconvolution with weighted robust linear regression',
                                   'Regression-based deconvolution with robust partial correlations'))
  }else if(show_description ==F){
    l = data.frame(method = c('nnls','cibersort','MuSiC','wRLM','RPC'),
                   function_to_call = c('deconv_regression_nnls',
                                        'deconv_regression_cibersort',
                                        'deconv_regression_MuSiC',
                                        'deconv_regression_wRLM',
                                        'deconv_RPC'),
                   suggested_packages = c('nnls','e1071/parallel/preprocessCore','MuSiC/SingleCellExperiment','LinDeconSeq','EpiDISH'))
  }
  return(l)
}

#' List available methods to perform marker-based deconvolution
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#'
#' @return a dataframe of available marker-based deconvolution methods along with the corresponding function names to execute these methods.
#' @export
#'
#' @examples
#' list_deconv_marker()
list_deconv_marker = function(show_description = F){
  if(show_description == T){
    l = data.frame(method = c('firstPC','gsva','debCAM','TOAST'),
                   function_to_call = c('deconv_marker_firstPC', 'deconv_marker_gsva', 'deconv_marker_debCAM', 'deconv_marker_TOAST'),
                   suggested_packages = c('FactoMineR','GSVA','debCAM','TOAST'),
                   description = c('Marker-based deconvolution with first principal score (abundance estimate only)',
                                   'Marker-based deconvolution with GSVA (abundance estimate only)',
                                   'Marker-based deconvolution with debCAM::AfromMarkers()',
                                   'Marker-based deconvolution with TOAST::MDeconv()'))
  }else if(show_description ==F){
    l = data.frame(method = c('firstPC','gsva','debCAM','TOAST'),
                   function_to_call = c('deconv_marker_firstPC', 'deconv_marker_gsva', 'deconv_marker_debCAM', 'deconv_marker_TOAST'),
                   suggested_packages = c('FactoMineR','GSVA','debCAM','TOAST'))
  }
  return(l)
}


#' List available methods to perform reference-free deconvolution
#'
#' @param show_description a logical variable to determine whether to show full description of each method
#'
#' @return a dataframe of available reference-free deconvolution methods along with the corresponding function names to execute these methods.
#' @export
#'
#' @examples
#' list_deconv_refFree()
list_deconv_refFree = function(show_description = F){
  if(show_description == T){
    l = data.frame(method = c('linseed','debCAM'),
                   function_to_call = c('deconv_refFree_linseed', 'deconv_refFree_debCAM'),
                   suggested_packages = c('linseed','debCAM'),
                   description = c('Reference free deconvolution with linseed',
                                   'Reference free deconvolution with debCAM'))
  }else if(show_description ==F){
    l = data.frame(method = c('linseed','debCAM'),
                   function_to_call = c('deconv_refFree_linseed', 'deconv_refFree_debCAM'),
                   suggested_packages = c('linseed','debCAM'))
  }
  return(l)
}




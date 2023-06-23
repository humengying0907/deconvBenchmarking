
#' Bayesian-based deconvolution with InstaPrism
#'
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param scExpr scRNA expression to build the reference, with genes in rows and samples in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_cellState column name that corresponds to cellState in scMeta. Set to NA if not available
#' @param colnames_of_sample 	column name that corresponds to sampleID in scMeta
#' @param key InstaPrism argument: name of the malignant cell type. Upon setting the key parameter, the updated malignant reference
#'    will be unique for each individual. Set to NA if there is no malignant cells in the problem, and the updated reference will be the same for all the individuals
#' @param n.iter InstaPrism argument: number of iterations. Default = 100
#' @param n.core 	number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results using both initial reference and updated reference.  If the 'key' argument is valid,
#'    it also includes the results using the malignant_updated reference (with malignant updated only)
#' @export
deconv_Bayesian_InstaPrism = function(bulk_expr,
                                      scExpr, scMeta,
                                      colnames_of_cellType = NA,
                                      colnames_of_cellState = NA,
                                      colnames_of_sample = NA,
                                      key = NA,
                                      n.iter = 100,
                                      n.core = 1){
  require(InstaPrism,quietly = T)

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cell_type in scMeta')
  }

  cell_type_labels = scMeta[,colnames_of_cellType]

  if(is.na(colnames_of_cellState)){
    if(is.na(colnames_of_sample)){
      cell_state_labels = scMeta[,colnames_of_cellType]
    }else{

      if(!is.na(key)){

        key_id = which(scMeta[,colnames_of_cellType] == key)
        cell_state_labels = scMeta[,colnames_of_cellType]
        cell_state_labels[key_id] = paste0(key,'_',scMeta[,colnames_of_sample][key_id])
      }

    }
  }else{
    cell_state_labels = scMeta[,colnames_of_cellState]
  }

  refPhi_obj = refPrepare(scExpr, cell_type_labels, cell_state_labels)

  InstaPrism_res = list()
  # InstaPrism with iniitial phi
  InstaPrism.res.initial = InstaPrism(input_type = 'refPhi',bulk_Expr = bulk_expr,refPhi = refPhi_obj, n.iter = n.iter,  n.core =  n.core)
  InstaPrism_res$Bayesian_InstaPrism_initial = InstaPrism.res.initial@Post.ini.ct@theta

  # InstaPrism with malignant reference updated only
  if(!is.na(key)){

    updatedMal_obj = InstaPrism_update(InstaPrism.res.initial,
                                        bulk_Expr = bulk_expr,
                                        n.iter = n.iter,
                                        cell.types.to.update = NULL,
                                        key = key,
                                        keep.phi = 'phi.ct')

    InstaPrism_res$Bayesian_InstaPrism_updatedMal = updatedMal_obj@theta
  }

  # InstaPrism with all cell types updated
  if(is.na(key)){
    cell.types.to.update = unique(cell_type_labels)
  }else{
    cell.types.to.update = unique(cell_type_labels)
    cell.types.to.update = cell.types.to.update[cell.types.to.update!=key]
  }

  updatedAll_obj  = InstaPrism_update(InstaPrism.res.initial,
                                      bulk_Expr = bulk_expr,
                                      n.iter = n.iter,
                                      cell.types.to.update = cell.types.to.update,
                                      key = key,
                                      keep.phi = 'phi.ct')

  InstaPrism_res$Bayesian_InstaPrism_updatedAll = updatedAll_obj@theta

  return(InstaPrism_res)
}


#' Bayesian-based deconvolution
#'
#' @param methods a character vector indicating which methods to use. Use list_deconv_Bayesian() to check for available method names.
#' @param bulk_expr bulk expression to deconvolute, with genes in rows and samples in columns
#' @param scExpr scRNA expression to build the reference, with genes in rows and samples in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_cellState column name that corresponds to cellState in scMeta. Set to NA if not available
#' @param colnames_of_sample 	column name that corresponds to sampleID in scMeta
#' @param key InstaPrism argument: name of the malignant cell type. Upon setting the key parameter, the updated malignant reference
#'    will be unique for each individual. Set to NA if there is no malignant cells in the problem, and the updated reference will be the same for all the individuals
#' @param n.iter InstaPrism argument: number of iterations. Default = 100
#' @param n.core 	number of cores to use for parallel programming. Default = 1
#'
#' @return a list of Bayesian-based deconvolution results
#' @export
#'
#' @examples
#' \dontrun{
#' deconv_Bayesian(methods = c('InstaPrism'),
#'                 bulk_expr = bulk_expr,
#'                 scExpr = scExpr,
#'                 scMeta = scMeta,
#'                 colnames_of_cellType = 'cell_type',
#'                 colnames_of_sample = 'sampleID',
#'                 key = 'malignant',
#'                 n.core = 4)
#'
#' }
deconv_Bayesian = function(methods,
                           bulk_expr,
                           scExpr, scMeta,
                           colnames_of_cellType = NA,
                           colnames_of_cellState = NA,
                           colnames_of_sample = NA,
                           key = NA,
                           n.iter = 100,
                           n.core = 1){

  l = list_deconv_Bayesian(show_description = T)

  Bayesian_deconvRes = list()

  for(method in methods){
    description = l$description[l$method == method]
    message(description)

    switch(method,
           InstaPrism = {
             res = deconv_Bayesian_InstaPrism(bulk_expr,
                                              scExpr, scMeta,
                                              colnames_of_cellType,
                                              colnames_of_cellState,
                                              colnames_of_sample,
                                              key,
                                              n.iter,
                                              n.core)
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_deconv_marker() to check for available methods"))
             result <- NULL
           })

    Bayesian_deconvRes = c(Bayesian_deconvRes,res)

  }
  return(Bayesian_deconvRes)
}





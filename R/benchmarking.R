
#' Training/testing cells splitting
#'
#' @param scMeta single cell annotation file with cells as rownames
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param training_ratio a numeric value indicating ratio of training cells. Default = 0.5
#' @param split_by a string vector indicating how to split training and testing cells. Available options includes 'cell' and 'sample'.
#'    When split_by = 'cell', for each cell type, a pre-defined ratio of cells will be sampled as training cells. When split = 'sample',
#'    a pre-defined ratio of samples will be sampled as training samples, and all the cells associated will become training cells.
#'
#' @return a list of training cells and testing cells
#' @export
#'
#' @examples
#'\dontrun{
#'# when splitting by cells
#'create_train_test_splitting(scMeta = scMeta, colnames_of_cellType = 'cell_type',
#'                            training_ratio = 0.5, split_by = 'cell')
#'
#'# when splitting by samples
#'create_train_test_splitting(scMeta = scMeta, colnames_of_sample = 'sampleID',
#'                            training_ratio = 0.5, split_by = 'sample')
#'}
create_train_test_splitting<-function(scMeta,
                                      colnames_of_cellType = NA,
                                      colnames_of_sample = NA,
                                      training_ratio = 0.5,
                                      split_by = 'sample'){

  if(split_by =='cell'){
    if(is.na(colnames_of_cellType)){
      stop('please provide "colnames_of_cellType" argument when split_by = "cell"')
    }

    scMeta_summary = scMeta %>% group_by_(colnames_of_cellType) %>%
      summarise(n=n()) %>% as.data.frame()

    scMeta_summary$training_n = ceiling(scMeta_summary$n * training_ratio)
    cell_types = unique(scMeta[,colnames_of_cellType])

    training_cells = c()
    for(i in 1:nrow(scMeta_summary)){
      cell_names = rownames(scMeta)[scMeta[,colnames_of_cellType] == scMeta_summary[i,1]]
      training_cells = c(training_cells,sample(cell_names,size = scMeta_summary[i,'training_n'],replace = F))
    }

    testing_cells = rownames(scMeta)[!rownames(scMeta) %in% training_cells]

    splitting_list=list()
    splitting_list[['training_cells']]=training_cells
    splitting_list[['testing_cells']]=testing_cells
  }else if(split_by == 'sample'){
    if(is.na(colnames_of_sample)){
      stop('please provide "colnames_of_sample" argument when split_by = "sample"')
    }

    unique_samples = unique(scMeta[,colnames_of_sample])
    training_samples = sample(unique_samples, size= ceiling(length(unique_samples)*training_ratio),replace = F)
    training_cells = rownames(scMeta)[scMeta[,colnames_of_sample] %in% training_samples]
    testing_cells = rownames(scMeta)[!rownames(scMeta) %in% training_cells]

    splitting_list=list()
    splitting_list[['training_cells']]=training_cells
    splitting_list[['testing_cells']]=testing_cells
  }else{
    stop('Invalid split_by argument provided')
  }
  return(splitting_list)
}

#' Create an object for deconvolution benchmarking
#' @description This function takes scRNA profile as input and generate an object intended for future deconvolution benchmarking.
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta.
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. This argument is required for 'semi', 'heter', and 'SCDC' bulk simulation methods;
#'    this argument is also required when split_by = 'sample'
#'
#' @param nbulk number of simulated bulk samples. This argument is required when simulated_frac is not available; in this case, the function
#'    will generate simulated_frac automatically using fracSimulator_Beta() with the nbulk argument.
#' @param fixed_cell_type argument for fracSimulator_Beta() function,  this argument will be called required when simulated_frac is not available.
#'    It denotes a character denotes the target cell type for which we strive to faithfully preserve its distribution. It is recommended to set this parameter
#'    to the name of the malignant cell types. If left undefined, the fracSimulator_Beta() will automatically select the most abundant cell type as 'fixed_cell_type'.
#'
#' @param training_ratio ratio of training cells. Default = 0.5
#' @param split_by a string vector indicating how to split training and testing cells. Available options includes 'cell' and 'sample'.
#'    When split_by = 'cell', for each cell type, a pre-defined ratio of cells will be sampled as training cells. When split = 'sample',
#'    a pre-defined ratio of samples will be sampled as training samples, and all the cells associated will become training cells.
#'
#' @param bulkSimulator_methods a character vector indicating which bulk simulation methods to use. Use list_bulkSimulator() to check for available method names.
#'    Set to NULL if no bulk simulation is needed.
#' @param colnames_of_subcluster column name that corresponds to subcluster info in scMeta, where subcluster contains sub-clustering information for each cellType.
#'    This is an argument required for 'heter_sampleIDfree' method only. Set to NA if subclustering information is not available; the function will generate subclustering information automatically
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns. If set to NULL, the function will generate simulated_frac automatically
#'    using 'nbulk' and 'fixed_cell_type' arguments.
#' @param heter_cell_type name of the cell_type to maintain the highest level of heterogeneity. It is recommended to set this parameter to the name of the malignant cell types.
#'     This argument is required for 'semi' and 'heter_sampleIDfree' bulk simulation methods
#' @param dirichlet_cs_par a numeric value determine the dispersion level of subclusters. With lower value indicating higher dispersion level. Default = 0.1. This is an argument required for 'heter_sampleIDfree' simulation method only.
#' @param min.subcluster.size minimum size of a subcluster. This argument is required when colnames_of_subcluster is NA. Default = 20. This is an argument for 'heter_sampleIDfree' simulation method only.
#'
#' @param refMarkers_methods a character vector specifying the desired methods for generating cell-type specific markers. Use list_refMarkers() to check for available method names.
#' @param colnames_of_cellState column name that corresponds to the cellState in scMeta. This argument is only required for scran-based marker identification,
#'    where the differential expression (DE) analysis is performed between every pair of cell states from different cell types. Set this parameter to NA if the information is not available.
#'    In that case, the function will treat all cells from the same cell type identically.
#'
#' @param refMatrix_methods a character vector specifying the desired methods for generating signature matrices. Use list_refMarix() to check for available method names.
#'
#' @param include_tcga a logicial variable determining whether to include tcga in the output object. If True, the function will download associated TCGA cohort from xnea browser
#' @param tcga_abbreviation a character indicating tcga abbreviation for the tcga cohort to include, for example 'SKCM'
#' @param purity_methods a character vector indicating tumor purity estimation method that is utilized as a means of estimating the malignant proportion within the exported object for TCGA expression data.
#'    Available methods include 'ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC' and 'CPE', Default = c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE')
#'
#' @param create_autogeneS_input a logical variable determine whether to create input data for autogeneS, which is a python based approach to construct signature matrix. If true, the function
#'    will automatically export the input data in 'autogeneS_input' folder and generate a command to run autogeneS with default or user-defined autogeneS hyperparameters
#' @param autogeneS_input_file_name desired file name to save input data for autogeneS
#' @param create_cibersortx_input a logical variable determine whether to create input data for cibersortx, which is a web-server to construct signature matrix. If true, the function will
#'    automatically export the input data in 'cibersortx_input' folder
#' @param cibersortx_input_file_name desired file name to save input data for cibersortx
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param ... additional arguments to be passed to the following functions: bulkSimulator(), refMarkers(), refMatrix(),
#'    pre_refMatrix_autogeneS() and pre_refMatrix_cibersortx()
#'
#'
#'
#' @details
#' This function takes scRNA profile as input and generate an object intended for future deconvolution benchmarking.
#'    The function performs the following steps: (1) It divides the cells into training and testing cells;
#'    (2) the training cells are utilized to generate reference profiles, such as markers and signature matrices;
#'    (3) the testing cells are used to generate simulated bulk expression, which is then employed for deconvolution purposes;
#'    (4) additionally, this function offers the flexibility to include a TCGA cohort as part of the object for future deconvolution benchmarking
#'
#' @section Training/testing splitting:
#' This function first splits the cells into training cells and testing cells using \code{\link{create_train_test_splitting}}. Important parameters
#' controlling how the dataset is splitted including:
#' \itemize{
#' \item \code{training_ratio}, which specifies the percentage of cells used as training cells
#' \item \code{split_by}, which determines whether to split the cells randomly or constrained on samples
#' }
#'
#' @section Reference profiles from training cells:
#' Training cells generated from splitting step will be utilized to generate reference profiles, such as markers and signature matrices. Specifically,
#' the exported \code{marker_list} in the function output is generated from \code{\link{refMarkers}}, and the exported \code{sigMarix_list} is
#' generated using \code{\link{refMatrix}}. Important parameters in reference construction including:
#' \itemize{
#' \item \code{refMarkers_methods}
#' \item \code{refMatrix_methods}
#' }
#' Set the above parameters to NULL if the associated reference profiles is not needed.
#'
#' @section Bulk simulation using testing cells:
#' Testing cells are used to generate simulated bulk expression using \code{\link{bulkSimulator}}. By default \code{bulkSimulator()} requires a pre-defined
#' fraction matrix \code{simulated_frac} determining cellular proportions to be aggregated. When \code{simulated_frac} is not available, the
#' function will generate a \code{simulated_frac} automatically using \code{\link{fracSimulator_Beta}}. Important parameters in bulk simulation including:
#' \itemize{
#' \item \code{bulkSimulator_methods}, which specifies which bulk simulation strategies to use
#' }
#'
#' @section Export files for autogeneS and cibersortx:
#' This function also provides the option to prepare input for autogeneS and cibersortx using the training cells. For details about the preparation steps and additional arguments associated,
#' check \code{\link{pre_refMatrix_autogeneS}} and \code{\link{pre_refMatrix_cibersortx}}.
#'
#' @return a list containing the following elements: 1) a list of training/testing cells; 2) a list of simulated bulk object and/or tcga expression;
#'    3) a list of cell-type specific markers; 4) a list of signature matrices
#' @export
#'
#' @examples
#' \dontrun{
#' # a standard benchmarking pipeline
#' benchmarking_init(scExpr = scExpr,
#'                   scMeta = scMeta,
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # argument for fracSimulator_Beta()
#'                   nbulk = 100,
#'                   fixed_cell_type = 'malignant',
#'
#'                   # argument for training/testing splitting:
#'                   training_ratio = 0.5,
#'                   split_by = "cell",
#'
#'                   # arguments for bulk simulation
#'                   bulkSimulator_methods = c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
#'
#'                   # argument for semi/heter_sampleIDfree bulk simulation method
#'                   heter_cell_type = 'malignant',
#'                   dirichlet_cs_par = 0.1,
#'                   min.subcluster.size = 20,
#'
#'                   # argument for marker constructions
#'                   refMarkers_methods = c('limma','scran'),
#'
#'                   # arguments for signature matrics construction
#'                   refMatrix_methods = c('raw','limma','scran'),
#'
#'                   # arguments to include tcga
#'                   include_tcga = T,
#'                   tcga_abbreviation = 'SKCM',
#'                   purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE'),
#'
#'                   # export files for autogeneS and cibersortx
#'                   create_autogeneS_input = T,
#'                   create_cibersortx_input = T
#'
#'                   n.core = 4
#'                   )
#'
#' # generate a benchmarking object containing only TCGA cohort
#' # and use all the single cells to generate reference markers and signature matrices
#' benchmarking_init(scExpr = scExpr,
#'                   scMeta = scMeta,
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # use all cells to build scRNA reference
#'                   training_ratio = 1,
#'
#'                   # set bulkSimulator_methods to NULL to disable bulk simulation
#'                   bulkSimulator_methods = NULL,
#'
#'                   # argument for marker constructions
#'                   refMarkers_methods = c('limma','scran'),
#'
#'                   # arguments for signature matrics construction
#'                   refMatrix_methods = c('raw','limma','scran'),
#'
#'                   # arguments to include tcga
#'                   include_tcga = T,
#'                   tcga_abbreviation = 'SKCM',
#'                   purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE')
#'                   )
#' }
benchmarking_init = function(scExpr,scMeta,

                             # arguments required in multiple steps
                             colnames_of_cellType = NA,
                             colnames_of_sample = NA,

                             # arguments for fracSimulator_Beta()
                             nbulk = 100,
                             fixed_cell_type = NA,

                             # arguments for create_train_test_splitting()
                             training_ratio = 0.5,
                             split_by = 'cell',

                             # arguments for bulkSimulator()
                             bulkSimulator_methods = NULL,
                             colnames_of_subcluster = NA,
                             simulated_frac = NULL,
                             heter_cell_type = NA,
                             dirichlet_cs_par = 0.1,
                             min.subcluster.size = 20,

                             # arguments for refMarkers
                             refMarkers_methods = c('limma','scran'),
                             colnames_of_cellState = NA,

                             # arguments for refMatrix
                             refMatrix_methods = c('raw','limma','scran'),

                             # arguments for including TCGA as part of evaluation
                             include_tcga = F,
                             tcga_abbreviation = NA,
                             purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE'),

                             # arguments for pre_refMatrix_autogeneS
                             create_autogeneS_input = F,
                             autogeneS_input_file_name = NULL,

                             # arguments for pre_refMatrix_cibersortx
                             create_cibersortx_input = F,
                             cibersortx_input_file_name = NULL,

                             n.core = 1,
                             ...
){

  arg_list = list(...)
  fracSimulator_Beta_args = arg_list[names(arg_list) %in% c('min.frac','showFractionPlot')]
  bulkSimulator_args = arg_list[names(arg_list) %in% c('ncells_perSample','min_chunkSize','use_chunk','min.percentage','max.percentage',
                                                       'seed','use_simulated_frac_as_prop_mat','disease','ct.sub','samplewithRep')]
  refMarkers_args = arg_list[names(arg_list) %in% c('hv_genes','log2FC','log2FC_flexible','minimum_n','maximum_n','spec_cutoff_for_DE','sigMatrixList')]
  refMatrix_args = arg_list[names(arg_list) %in% c('hv_genes','log2FC','log2FC_flexible','minimum_n','maximum_n','spec_cutoff_for_DE','markerList ')]
  pre_refMatrix_autogeneS_args = arg_list[names(arg_list) %in% c('hv_genes','max.spec_cutoff_for_autogeneS','display_autogeneS_command','ngen',
                                                                 'seed','nfeatures','mode')]
  pre_refMatrix_cibersortx_args = arg_list[names(arg_list) %in% c('downsample','downsample_ratio')]

  benchmarking_obj = list()

  message('split scRNA into training and testing')
  splitting = create_train_test_splitting(scMeta, colnames_of_cellType, colnames_of_sample,
                                          training_ratio, split_by)

  # in case that training_cells and testing_cells belong to non-identical cell_types,
  # which may occur when certain cell_types are only present in a limited number of samples

  scMeta_train = scMeta[splitting$training_cells,]
  scMeta_test = scMeta[splitting$testing_cells,]

  overlapped_celltypes = intersect(scMeta_train[,colnames_of_cellType],
                                   scMeta_test[,colnames_of_cellType])

  # reorganize scMeta and scExpr using overlapped_celltypes
  scMeta_renamed = data.frame(row.names = rownames(scMeta),
                              sampleID = scMeta[,colnames_of_sample],
                              cell_type = scMeta[,colnames_of_cellType])

  scMeta_renamed = scMeta_renamed[scMeta_renamed$cell_type %in% overlapped_celltypes,]

  splitting$training_cells = splitting$training_cells[splitting$training_cells %in% rownames(scMeta_renamed)]
  splitting$testing_cells = splitting$testing_cells[splitting$testing_cells %in% rownames(scMeta_renamed)]

  # make training and testing profiles
  scExpr_train = scExpr[,splitting$training_cells]
  scMeta_train = scMeta[splitting$training_cells,]

  scExpr_test = scExpr[,splitting$testing_cells]
  scMeta_test = scMeta[splitting$testing_cells,]

  if(is.na(colnames_of_cellType)){
    cell_type_labels_train = NULL
  }else{
    cell_type_labels_train = scMeta_train[,colnames_of_cellType]
  }

  if(is.na(colnames_of_cellState)){
    cell_state_labels_train = NULL
  }else{
    cell_state_labels_train = scMeta_train[,colnames_of_cellState]
  }

  if(!is.null(bulkSimulator_methods)){

    if(is.null(simulated_frac)){
      message('simulate realistic cell-type fractions')

      simulated_frac = do.call(fracSimulator_Beta,c(list(scMeta = scMeta_renamed,
                                                         n = nbulk,
                                                         colnames_of_sample  = 'sampleID',
                                                         colnames_of_cellType  = 'cell_type',
                                                         fixed_cell_type = fixed_cell_type),
                                                    fracSimulator_Beta_args))


      cat('\n')
    }

    message('generate a list of bulk simulation objects using testing scRNA')
    bulk_simulation_obj = do.call(bulkSimulator,c(list(methods = bulkSimulator_methods,
                                                       scExpr = scExpr_test,
                                                       scMeta = scMeta_test,
                                                       colnames_of_cellType = colnames_of_cellType,
                                                       colnames_of_subcluster = colnames_of_subcluster,
                                                       colnames_of_sample = colnames_of_sample,
                                                       simulated_frac = simulated_frac,
                                                       heter_cell_type  = heter_cell_type,
                                                       dirichlet_cs_par = dirichlet_cs_par ,
                                                       min.subcluster.size = min.subcluster.size,
                                                       nbulk = nbulk,
                                                       n.core = n.core),
                                                  bulkSimulator_args))


    cat('\n')

  }else{
    bulk_simulation_obj = list()
  }

  if(include_tcga == T){
    message('include TCGA in the bulk object for further evaluation')
    to_use_genes = rownames(scExpr)
    bulk_simulation_obj[[paste0('tcga_',tcga_abbreviation)]] = build_tcga_obj(tcga_abbreviation,
                                                                              to_use_genes,
                                                                              purity_methods)
    cat('\n')
  }

  if(!is.null(refMarkers_methods)){
    message('generate a list of cell-type specific markers from training scRNA')
    marker_list = do.call(refMarkers,c(list(methods = refMarkers_methods,
                                            scExpr = scExpr_train,
                                            cell_type_labels = cell_type_labels_train,
                                            cell_state_labels  = cell_state_labels_train),
                                       refMarkers_args))


    cat('\n')
  }else{
    marker_list = list()
  }

  if(!is.null(refMarkers_methods)){
    message('generate a list of signature matrices from training scRNA')
    query = all(intersect(refMarkers_methods,refMatrix_methods) %in% c('limma','scran')) & length(intersect(refMarkers_methods,refMatrix_methods)>0)

    if(query == T){
      refMatrix_methods_remaining = refMatrix_methods[!refMatrix_methods %in% refMarkers_methods]
      sigMarix_list = do.call(refMatrix,c(list(methods = c(refMatrix_methods_remaining,'markerList'),
                                               scExpr = scExpr_train,
                                               cell_type_labels = cell_type_labels_train,
                                               markerList = marker_list),
                                          refMatrix_args))

    }else{
      sigMarix_list = do.call(refMatrix,c(list(methods = refMatrix_methods,
                                               scExpr = scExpr_train,
                                               cell_type_labels = cell_type_labels_train,
                                               cell_state_labels  = cell_state_labels_train),
                                          refMatrix_args))

    }
    cat('\n')
  }else{
    sigMarix_list = list()
  }


  if(create_autogeneS_input==T){
    message('prepare input for autogeneS using training scRNA')
    do.call(pre_refMatrix_autogeneS,c(list(scExpr = scExpr_train,
                                           cell_type_labels = cell_type_labels_train,
                                           autogeneS_input_file_name = autogeneS_input_file_name),
                                      pre_refMatrix_autogeneS_args))


    cat('\n')
  }

  if(create_cibersortx_input ==T){
    message('prepare input for cibersortx using training scRNA')
    do.call(pre_refMatrix_cibersortx,c(list(scExpr = scExpr_train,
                                            cell_type_labels = cell_type_labels_train,
                                            scMeta = scMeta_train,
                                            colnames_of_sample = colnames_of_sample,
                                            colnames_of_cellType = colnames_of_cellType,
                                            cibersortx_input_file_name = cibersortx_input_file_name),
                                       pre_refMatrix_cibersortx_args))

        cat('\n')
  }

  # organize all together
  benchmarking_obj$splitting = splitting
  benchmarking_obj$bulk_simulation_obj = bulk_simulation_obj
  benchmarking_obj$marker_list = marker_list
  benchmarking_obj$sigMarix_list = sigMarix_list

  return(benchmarking_obj)
}



#' Perform deconvolution on the benchmarking_obj
#'
#' @param benchmarking_obj a benchmarking_obj generated from benchmarking_init() function
#'
#' @param marker_based_methods a character vector indicating which methods to use for marker_based deconvolution.
#'    Use list_deconv_marker() to check for available method names. Set to NULL if no marker_based deconvolution is needed.
#' @param alpha TOAST::MDeconv() argument. A vector including the prior mean for all cell types. Set to NULL to for TOAST/-P (partial reference-free deconvolution without prior)
#' @param sigma TOAST::MDeconv() argument. A vector including the prior standard deviation for all cell types.  Set to NULL to for TOAST/-P (partial reference-free deconvolution without prior)
#' @param epsilon TOAST::MDeconv() argument. A numeric variable to control the level of convergence. With a large epsilon, it converges quickly and the algorithm may not converge well.
#'    With a small epsilon, it converges slower and the algorithm may converge "too much". The default value is 1e-3, which we find is a reasonable threshold.
#' @param maxIter TOAST::MDeconv() argument. Number of maximum iterations.
#'
#' @param regression_based_methods a character vector indicating which methods to use for regression_based deconvolution.
#'    Use list_deconv_regression() to check for available method names. Set to NULL if no regression_based deconvolution is needed.
#'
#' @param cibersort_path path to CIBERSORT.R. This argument is required for 'cibersort' method
#' @param QN_cibersort a logical variable determining whether to quantile normalize the input matrix. Default = F. This argument is required for 'cibersort' method
#' @param skip_raw_cibersort a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources. This argument is required for 'cibersort' method
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns. This argument is required for 'MuSiC' method
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr. This argument is required for 'MuSiC' method
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta. This argument is required for 'MuSiC' method
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. This argument is required for 'MuSiC' method
#' @param normalize_MuSiC MuSiC::music_prop() normalize argument. Default = F. This argument is required for 'MuSiC' method
#'
#' @param skip_raw_wRLM a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources. This argument is required for 'wRLM' method
#' @param weight_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to weight gene or not. Default = T. This argument is required for 'wRLM' method
#' @param intercept_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to add intercept when using robust linear regression. Default = T. This argument is required for 'wRLM' method
#' @param scale_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to scale bulk gene expression or not. Default = T. This argument is required for 'wRLM' method
#' @param QN_wRLM LinDeconSeq::deconSeq() argument, a logical variable determining whether to quantile normalize bulk expression profile. Default = T. This argument is required for 'wRLM' method
#'
#' @param skip_raw_RPC a logical variable indicating whether to skip the signature matrix named 'raw', which corresponds to the raw reference matrix without gene filtering.
#'    Setting this argument to TRUE is recommended to save computation time and resources. This argument is required for 'RPC' method
#' @param maxit_RPC EpiDISH::epidish() argument, an integer indicating the limit of the number of IWLS iterations. Default = 100. This argument is required for RPC method
#'
#' @param refFree_methods a character vector indicating which methods to use for reference free deconvolution. Use list_deconv_refFree() to check for available method names.
#' @param k argument for reference free methods: number of cell types in bulk expression. If set to NA, will be default to true number of cell types in bulk expression
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
#' @param Bayesian_methods a character vector indicating which Bayesian_methods to use. Use list_deconv_Bayesian() to check for available method names.
#' @param colnames_of_cellState column name that corresponds to cellState in scMeta. Set to NA if not available
#' @param key InstaPrism argument: name of the malignant cell type. Upon setting the key parameter, the updated malignant reference
#'    will be unique for each individual. Set to NA if there is no malignant cells in the problem, and the updated reference will be the same for all the individuals
#' @param n.iter InstaPrism argument: number of iterations. Default = 100
#'
#' @param immunedeconv_methods a character vector indicating which methods to use from immunedeconv package. Available methods include:
#'    'xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'
#' @param tcga_abbreviation a character string indicating tcga-abbreviation for the bulk data to be deconvoluted, for example 'skcm'. Required for 'timer' and 'consensus_tme'.
#' @param tumor a logical variable to determine whether to use a signature matrix/procedure optimized for tumor samples, if supported by the method. Currently affects EPIC and quanTIseq.
#' @param scale_mrna logical. If FALSE, disable correction for mRNA content of different cell types. This is supported by methods that compute an absolute score (EPIC and quanTIseq)
#'
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list of deconvolution results for each bulk expression in the benchmarking_obj
#' @export
#'
#' @examples
#' \dontrun{
#' # first create a benchmarking_obj using benchmarking_init():
#' benchmarking_obj = benchmarking_init(scExpr = scExpr,
#'                   scMeta = scMeta,
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # argument for fracSimulator_Beta()
#'                   nbulk = 100,
#'                   fixed_cell_type = 'malignant',
#'
#'                   # argument for training/testing splitting:
#'                   training_ratio = 0.5,
#'                   split_by = "cell",
#'
#'                   # arguments for bulk simulation
#'                   bulkSimulator_methods = c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
#'
#'                   # argument for semi/heter_sampleIDfree bulk simulation method
#'                   heter_cell_type = 'malignant',
#'                   dirichlet_cs_par = 0.1,
#'                   min.subcluster.size = 20,
#'
#'                   # argument for marker constructions
#'                   refMarkers_methods = c('limma','scran'),
#'
#'                   # arguments for signature matrics construction
#'                   refMatrix_methods = c('raw','limma','scran'),
#'
#'                   # arguments to include tcga
#'                   include_tcga = T,
#'                   tcga_abbreviation = 'SKCM',
#'                   purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE'),
#'
#'                   # export files for autogeneS and cibersortx
#'                   create_autogeneS_input = T,
#'                   create_cibersortx_input = T
#'
#'                   n.core = 4
#'                   )
#'
#' # perform deconvolution on the benchmarking_obj using all categories of deconvolution methods, which includes:
#' # 1) marker-based methods: 'firstPC','gsva','debCAM','TOAST'
#' # 2) regression-based methods: 'nnls','cibersort','MuSiC','wRLM','RPC'
#' # 3) reference free methods: 'linseed','debCAM'
#' # 4) Bayesian-based methods: 'InstaPrism
#' # 5) other deconvolution methods from immunedeconv package:
#' #     'xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'
#' benchmarking_deconv(benchmarking_obj,
#'
#'                     # arguments for marker based methods
#'                     marker_based_methods = c('firstPC','gsva','debCAM','TOAST'),
#'
#'                     # arguments for regression-based methods
#'                     regression_based_methods =  c('nnls','cibersort','MuSiC','wRLM','RPC'),
#'                     cibersort_path = 'scripts/', # argument required for 'cibersort' method
#'                     scExpr = scExpr, # arguments required for 'MuSiC'
#'                     scMeta = scMeta,
#'                     colnames_of_cellType = 'cell_type',
#'                     colnames_of_sample = 'sampleID',
#'
#'                     # arguments for reference-free methods
#'                     refFree_methods = c('linseed','debCAM'),
#'
#'                     # arguments for Bayesian-based methods
#'                     Bayesian_methods = c('InstaPrism'),
#'                     key = 'malignant',  # this argument togther with 'colnames_of_sample' is highly recommended to run Bayesian based methods
#'
#'                     # arguments for other deconvolution methods
#'                     immunedeconv_methods = c('xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'),
#'                     tcga_abbreviation = 'SKCM', # arguments required for 'timer' and 'consensus_tme'
#'
#'
#'                     n.core = 4)
#'
#' # if only marker-based deconvolution and regression-based deconvolution is needed
#' benchmarking_deconv(benchmarking_obj,
#'                     marker_based_methods = c('firstPC','gsva','debCAM','TOAST'),
#'
#'                     regression_based_methods =  c('nnls','cibersort','MuSiC','wRLM','RPC'),
#'                     cibersort_path = 'scripts/',  # argument required for 'cibersort' method
#'                     scExpr = scExpr, # arguments required for 'MuSiC'
#'                     scMeta = scMeta,
#'                     colnames_of_cellType = 'cell_type',
#'                     colnames_of_sample = 'sampleID',
#'
#'                     # set methods argument to NULL to disable deconvolution with these methods
#'                     refFree_methods = NULL,
#'                     Bayesian_methods = NULL,
#'                     immunedeconv_methods = NULL,
#'
#'                     n.core = 4)
#'
#' }
benchmarking_deconv = function(benchmarking_obj,
                               # arguments for deconv_marker()
                               marker_based_methods = c('firstPC','gsva','debCAM','TOAST'),
                               alpha = NULL, sigma = NULL,
                               epsilon = 0.001, maxIter = 1000,

                               # arguments for deconv_regression()
                               regression_based_methods =  c('nnls','cibersort','MuSiC','wRLM','RPC'),
                               cibersort_path = NULL,
                               QN_cibersort = F, skip_raw_cibersort = T,
                               scExpr = NULL, scMeta = NULL, colnames_of_cellType = NA, colnames_of_sample = NA, normalize_MuSiC = F,
                               skip_raw_wRLM = TRUE, weight_wRLM = TRUE, intercept_wRLM = TRUE, scale_wRLM = FALSE, QN_wRLM = FALSE,
                               skip_raw_RPC = TRUE, maxit_RPC = 100,

                               # arguments for deconv_refFree()
                               refFree_methods = c('linseed','debCAM'),
                               k = NA,
                               corner.strategy = 2, dim.rdc = 10,
                               thres.low = 0.05, thres.high = 0.95, cluster.method = "K-Means", cluster.num = 50, MG.num.thres = 20,
                               lof.thres = 0.02, quickhull = TRUE, quick.select = NULL,
                               sample.weight = NULL, appro3 = TRUE, generalNMF = FALSE,
                               cores = NULL,

                               # arguments for deconv_Bayesian
                               Bayesian_methods = c('InstaPrism'),
                               colnames_of_cellState = NA,
                               key = NA,
                               n.iter = 100,

                               # arguments for methods from immunedeconv package
                               immunedeconv_methods = c('xcell','mcp_counter','epic','quantiseq','timer','abis',
                                                 'consensus_tme','estimate'),
                               tcga_abbreviation = NA,
                               tumor = T,
                               scale_mrna = T,

                               n.core = 1

){
  deconvResults = list()

  if('cibersort' %in% regression_based_methods){
    if(is.null(cibersort_path)){
      stop('please provide path to cibersort R script to run cibersort')
    }
  }

  if('MuSiC' %in% regression_based_methods){
    if(is.null(scExpr) | is.null(scMeta) | is.na(colnames_of_cellType) | is.na(colnames_of_sample)){
      stop('please provide the following required arguments to run MuSiC: scExpr, scMeta, colnames_of_cellType, colnames_of_sample')
    }
  }

  if(!is.null(immunedeconv_methods)){
    require(immunedeconv)
  }

  if('consensus_tme' %in% immunedeconv_methods | 'timer' %in% immunedeconv_methods){
    if(is.na(tcga_abbreviation)){
      stop('please provide tcga_abbreviation to run consensus_tme and/or timer')
    }
  }

  if('InstaPrism' %in% Bayesian_methods){
    if(is.null(scExpr) | is.null(scMeta) | is.na(colnames_of_cellType) ){
      stop('please provide the following required arguments to run InstaPrism: scExpr, scMeta, colnames_of_cellType;
           additionally, it is highly recommended to provide the "colnames_of_sample" and "key" arguments to run InstaPrism')
    }
  }

  for(bulk in names(benchmarking_obj$bulk_simulation_obj)){

    bulk_expr = benchmarking_obj[['bulk_simulation_obj']][[bulk]][['simulated_bulk']]
    bulk_deconvRes = list()

    if(!is.null(marker_based_methods)){
      marker_based_deconvRes = deconv_marker(marker_based_methods,
                                             bulk_expr,
                                             benchmarking_obj$marker_list,
                                             alpha, sigma,
                                             epsilon, maxIter,
                                             n.core)
      bulk_deconvRes = c(bulk_deconvRes,marker_based_deconvRes)
    }

    if(!is.null(regression_based_methods)){

      if('MuSiC' %in% regression_based_methods){
        scExpr_train = scExpr[,benchmarking_obj$splitting$training_cells]
        scMeta_train = scMeta[benchmarking_obj$splitting$training_cells,]
      }

      regression_based_deconvRes = deconv_regression(regression_based_methods,
                                                     bulk_expr,
                                                     benchmarking_obj$sigMarix_list,
                                                     cibersort_path,
                                                     QN_cibersort, skip_raw_cibersort,
                                                     scExpr_train,scMeta_train,
                                                     colnames_of_cellType, colnames_of_sample, normalize_MuSiC,
                                                     skip_raw_wRLM, weight_wRLM , intercept_wRLM, scale_wRLM, QN_wRLM,
                                                     skip_raw_RPC, maxit_RPC,
                                                     n.core)
      bulk_deconvRes = c(bulk_deconvRes,regression_based_deconvRes)
    }

    if(!is.null(refFree_methods)){
      if(is.na(k)){
        k = ncol(benchmarking_obj[['bulk_simulation_obj']][[bulk]][['simulated_frac']])
      }

      refFree_deconvRes = deconv_refFree(refFree_methods, bulk_expr,k,
                                         corner.strategy = 2, dim.rdc = 10,
                                         thres.low = 0.05, thres.high = 0.95, cluster.method = "K-Means", cluster.num = 50, MG.num.thres = 20,
                                         lof.thres = 0.02, quickhull = TRUE, quick.select = NULL,
                                         sample.weight = NULL, appro3 = TRUE, generalNMF = FALSE,
                                         cores = NULL)
      bulk_deconvRes = c(bulk_deconvRes,refFree_deconvRes)
    }

    if(!is.null(Bayesian_methods)){
      scExpr_train = scExpr[,benchmarking_obj$splitting$training_cells]
      scMeta_train = scMeta[benchmarking_obj$splitting$training_cells,]

      Bayesian_deconvRes = deconv_Bayesian(Bayesian_methods, bulk_expr,
                                           scExpr_train,scMeta_train,
                                           colnames_of_cellType,
                                           colnames_of_cellState,
                                           colnames_of_sample,
                                           key,
                                           n.iter,
                                           n.core)
      bulk_deconvRes = c(bulk_deconvRes,Bayesian_deconvRes)

    }

    if(!is.null(immunedeconv_methods)){
      immunedeconv_deconvRes = list()

      for(method in immunedeconv_methods){
        if(method %in% c('consensus_tme','timer')){
          indications = rep(tcga_abbreviation,ncol(bulk_expr))
          indications = stringr::str_to_lower(indications)
        }else{
          indications = NULL
        }

        immunedeconv_deconvRes[[method]] = immunedeconv::deconvolute(bulk_expr,method,
                                                                     indications,
                                                                     tumor = tumor,
                                                                     scale_mrna = scale_mrna)  %>% column_to_rownames('cell_type') %>% as.matrix()
      }
      bulk_deconvRes = c(bulk_deconvRes, immunedeconv_deconvRes)
    }

    deconvResults[[bulk]] = list(deconvRes = lapply(bulk_deconvRes,t),
                                 true_frac = benchmarking_obj[['bulk_simulation_obj']][[bulk]][['simulated_frac']])
    message(paste('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> finish deconvolution for',bulk,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))
    cat('\n')
  }
  return(deconvResults)
}


#' Evaluate deconvolution performance
#' @description This function evaluates deconvolution performance on the deconvResults_obj returned from benchmarking_deconv() function.
#' @param deconvResults_obj a deconvResults_obj returned from deconvResults_obj() function
#'
#' @return a list of performance evaluation statistics including per-cell type correlation and RMSE values, specified for each simulated bulk expression
#' @export
#'
#' @examples
#' \dontrun{
#' # a complete pipeline to benchmarking deconvolution methods using bulk data simulated under various strategies
#'
#' # first create a benchmarking_obj using benchmarking_init():
#' benchmarking_obj = benchmarking_init(scExpr = scExpr,
#'                   scMeta = scMeta,
#'
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # argument for fracSimulator_Beta()
#'                   fixed_cell_type = 'malignant',
#'
#'                   # arguments for bulk simulation
#'                   bulkSimulator_methods = c('homo','semi','heter','heter_sampleIDfree','favilaco','immunedeconv','SCDC'),
#'
#'                   # argument for semi/heter_sampleIDfree bulk simulation method
#'                   heter_cell_type = 'malignant',
#'                   dirichlet_cs_par = 0.1,
#'                   min.subcluster.size = 20,
#'
#'                   # argument for marker constructions
#'                   refMarkers_methods = c('limma','scran'),
#'
#'                   # arguments for signature matrics construction
#'                   refMatrix_methods = c('raw','limma','scran'),
#'
#'                   # arguments to include tcga
#'                   include_tcga = T,
#'                   tcga_abbreviation = 'SKCM',
#'                   purity_methods =  c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE'),
#'
#'                   # export files for autogeneS and cibersortx
#'                   create_autogeneS_input = T,
#'                   create_cibersortx_input = T
#'
#'                   n.core = 4
#'                   )
#'
#' # perform deconvolution on the benchmarking_obj using all categories of deconvolution methods, which includes:
#' # 1) marker-based methods: 'firstPC','gsva','debCAM','TOAST'
#' # 2) regression-based methods: 'nnls','cibersort','MuSiC','wRLM','RPC'
#' # 3) reference free methods: 'linseed','debCAM'
#' # 4) Bayesian-based methods: 'InstaPrism
#' # 5) other deconvolution methods from immunedeconv package:
#' #     'xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'
#'
#' benchmarking_deconv(benchmarking_obj,
#'
#'                     # arguments for marker based methods
#'                     marker_based_methods = c('firstPC','gsva','debCAM','TOAST'),
#'
#'                     # arguments for regression-based methods
#'                     regression_based_methods =  c('nnls','cibersort','MuSiC','wRLM','RPC'),
#'                     cibersort_path = 'scripts/', # argument required for 'cibersort' method
#'                     scExpr = scExpr, # arguments required for 'MuSiC'
#'                     scMeta = scMeta,
#'                     colnames_of_cellType = 'cell_type',
#'                     colnames_of_sample = 'sampleID',
#'
#'                     # arguments for reference-free methods
#'                     refFree_methods = c('linseed','debCAM'),
#'
#'                     # arguments for Bayesian-based methods
#'                     Bayesian_methods = c('InstaPrism'),
#'                     key = 'malignant',  # this argument togther with 'colnames_of_sample' is highly recommended to run Bayesian based methods
#'
#'                     # arguments for other deconvolution methods
#'                     immunedeconv_methods = c('xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'),
#'                     tcga_abbreviation = 'SKCM', # arguments required for 'timer' and 'consensus_tme'
#'
#'
#'                     n.core = 4)
#'
#' # evaluate the performance of each deconvolution methods:
#' benchmarking_evalu(deconvResults_obj)
#' }
benchmarking_evalu = function(deconvResults_obj){

  get_maxCor = function(x){
    return(x$summ$cor)
  }

  get_RMSE = function(x){
    return(x$summ$RMSE)
  }

  deconvPerformance = list()

  for(bulk in names(deconvResults_obj)){
    Y = deconvResults_obj[[bulk]]$true_frac
    Y = Y[,order(colnames(Y)),drop = F]

    bulk_performance = list()

    # detailed per-method evaluation statistics
    deconvEvalu = list()
    for(method in names(deconvResults_obj[[bulk]]$deconvRes)){

      E = deconvResults_obj[[bulk]]$deconvRes[[method]]
      E = E[,order(colnames(E)),drop = F]

      if(ncol(Y) != ncol(E)){
        # find max-correlated column in deconRes for each cell type
        maxCorName = c()
        for(ct in colnames(Y)){
          id = apply(cor(E,Y[,ct],use = 'pairwise.complete.obs'),2,which.max)
          maxCorName = c(maxCorName, colnames(E)[id])
        }

        E = E[,maxCorName,drop = F]
        colnames(E) = make.unique(colnames(E))
      }else{
        if(all.equal(colnames(Y),colnames(E)) == F){
          # find max-correlated column in deconRes for each cell type
          maxCorName = c()
          for(ct in colnames(Y)){
            id = apply(cor(E,Y[,ct],use = 'pairwise.complete.obs'),2,which.max)
            maxCorName = c(maxCorName, colnames(E)[id])
          }

          E = E[,maxCorName, drop = F]
          colnames(E) = make.unique(colnames(E))
        }else{
          E = E
        }
      }

      m1 = gather(Y %>% as.data.frame(),cell_type,true_frac)
      m2 = gather(E %>% as.data.frame(),maxCorName,estimate)
      M = cbind(m1,m2)

      # add summary statistics
      summ <- M %>%
        group_by(cell_type) %>%
        summarise(
          RMSE = caret::RMSE(true_frac, estimate,na.rm = T),
          cor = cor(true_frac,estimate,method = 'pearson',use = 'pairwise.complete.obs')) %>%
        mutate_if(is.numeric, round, digits=2) %>% as.data.frame() %>% column_to_rownames('cell_type')

      summ$maxCorName = M$maxCorName[match(rownames(summ),M$cell_type)]
      deconvEvalu[[method]] = list(M = M, summ = summ)
    }


    bulk_performance$maxCor = do.call(cbind, lapply(deconvEvalu,get_maxCor))
    rownames(bulk_performance$maxCor) = colnames(Y)

    bulk_performance$RMSE = do.call(cbind, lapply(deconvEvalu,get_RMSE))
    rownames(bulk_performance$RMSE) = colnames(Y)

    bulk_performance$deconvEvalu = deconvEvalu


    deconvPerformance[[bulk]] = bulk_performance
  }

  return(deconvPerformance)
}




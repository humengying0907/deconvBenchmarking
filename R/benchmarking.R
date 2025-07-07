
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
#' @param fixed_cell_type argument for fracSimulator_Beta() function,  this argument is applicable when simulated_frac is not available.
#'    It is a character denoting the target cell type for which we strive to faithfully preserve its distribution. It is recommended to set this parameter
#'    to the name of the malignant cell types. If left undefined, the fracSimulator_Beta() will automatically select the most abundant cell type as 'fixed_cell_type'.
#'
#' @param training_ratio ratio of training cells. Default = 0.5
#' @param split_by a string vector indicating how to split training and testing cells. Available options includes 'cell' and 'sample'.
#'    When split_by = 'cell', for each cell type, a pre-defined ratio of cells will be sampled as training cells. When split = 'sample',
#'    a pre-defined ratio of samples will be sampled as training samples, and all the cells belonging to training samples will become training cells.
#'
#' @param bulkSimulator_methods a character vector indicating which bulk simulation methods to use. Use \code{\link{list_bulkSimulator}} to check for available method names and suggested packages associated them.
#'    Make sure you have the required packages installed to use these methods.Set to NULL if no bulk simulation is needed.
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns. If set to NULL, the function will generate simulated_frac automatically
#'    using 'nbulk' and 'fixed_cell_type' arguments.
#' @param heter_cell_type name of the cell_type to maintain the highest level of heterogeneity. It is recommended to set this parameter to the name of the malignant cell-type.
#'     This argument is required for 'semi' and 'heter_sampleIDfree' bulk simulation methods
#' @param ncells_perSample number of cells to aggregate for each simulated bulk sample. This is an argument required for 'homo', 'semi', 'immunedeconv' and 'SCDC' methods
#' @param min_chunkSize minimum number of cells required to construct a particular cell-type component in the simulated bulk, such as requiring at least 20 cells for B cells, at least 20 cells for T cells, and so forth. This is an argument required for 'semi' and 'heter' methods
#' @param use_chunk a character indicating which cells to pool together for the a particular cell-type component. Default='all' other options include 'random'.
#'    When use_chunk = 'all', use all the cells belonging to the same patient for a given cell type to generate the certain cell type component in the simulated bulk;
#'    when use_chunk = 'random', randomly select 50-100% of the cells belonging to the same patient for a given cell type. This is an argument required for 'semi' and 'heter' methods
#' @param colnames_of_subcluster column name that corresponds to subcluster info in scMeta, where subcluster contains sub-clustering information for each cellType.
#'    This is an argument required for 'heter_sampleIDfree' method only. Set to NA if subclustering information is not available; the function will then generate subclustering information automatically using the 'min.subcluster.size' argument
#' @param export_cellUsage a logical variable determining whether to export cell names used to generate the simulated bulk. Default = F. This is an argument is only applicable to 'homo', 'semi', 'heter' and 'heter_sampleIDfree' methods
#'
#' @param refMarkers_methods a character vector specifying the desired methods for generating cell-type specific markers. Use \code{\link{list_refMarkers}} to check for available method names and suggested packages associated them.
#'    Make sure you have the required packages installed to use these methods. Set to NULL if cell-type specific markers are not required
#' @param colnames_of_cellState column name that corresponds to the cellState in scMeta. This argument is only required for scran-based marker identification,
#'    where the differential expression (DE) analysis is performed between every pair of cell states from different cell types. Set this parameter to NA if the information is not available.
#'    In that case, the function will treat all cells from the same cell type identically.
#'
#' @param refMatrix_methods a character vector specifying the desired methods for generating signature matrices. Use \code{\link{list_refMarix}} to check for available method names and suggested packages associated them.
#'    Make sure you have the required packages installed to use these methods. Set to NULL if signature matrices are not needed.
#'
#' @param include_tcga a logical variable determining whether to include tcga in the output object. If True, the function will download the specified TCGA cohort from the xena browser
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
#'
#' @return a list containing the following elements: 1) a list of training/testing cells; 2) a list of simulated bulk object;
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
#'                   # arguments for bulk simulation: select bulk simulation methods
#'                   bulkSimulator_methods = c('homo','semi','heter','heter_sampleIDfree'),
#'
#'                   # argument required for semi/heter_sampleIDfree bulk simulation method
#'                   heter_cell_type = 'malignant',
#'
#'                   # general simulation parameters
#'                   ncells_perSample = 500,
#'                   min_chunkSize = 20,
#'                   use_chunk = 'random',
#'
#'                   # parameters for adjusting heterogeneity in sample ID-free bulk simulation
#'                   dirichlet_cs_par = 0.1,
#'                   min.subcluster.size = 20,
#'                   max.num.cs = NA,
#'
#'                   # argument for marker constructions
#'                   refMarkers_methods = c('limma','scran'),
#'
#'                   # arguments for signature matrics construction
#'                   refMatrix_methods = c('raw','limma','scran'),
#'
#'                   # export files for autogeneS and cibersortx
#'                   create_autogeneS_input = T,
#'                   create_cibersortx_input = T,
#'
#'                   n.core = 4
#'                   )
#'
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
                             simulated_frac = NULL,
                             heter_cell_type = NA,
                             ncells_perSample=500,
                             min_chunkSize=5,
                             use_chunk = 'all',
                             colnames_of_subcluster = NA,
                             export_cellUsage = F,

                             # arguments for refMarkers
                             refMarkers_methods = c('limma','scran'),
                             colnames_of_cellState = NA,

                             # arguments for refMatrix
                             refMatrix_methods = c('raw','limma','scran'),

                             # arguments for including TCGA as part of evaluation
                             include_tcga = F,
                             tcga_abbreviation = NA,

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
  bulkSimulator_args = arg_list[names(arg_list) %in% c('dirichlet_cs_par','min.subcluster.size','max.num.cs',
                                                       'min.percentage','max.percentage','seed',
                                                       'disease','ct.sub','prop_mat','samplewithRep')]
  refMarkers_args = arg_list[names(arg_list) %in% c('hv_genes','log2FC','minimum_n','maximum_n','spec_cutoff_for_DE','sigMatrixList')]
  refMatrix_args = arg_list[names(arg_list) %in% c('hv_genes','log2FC','minimum_n','maximum_n','order_by','spec_cutoff_for_DE','markerList')]
  pre_refMatrix_autogeneS_args = arg_list[names(arg_list) %in% c('hv_genes','max.spec_cutoff_for_autogeneS','display_autogeneS_command','ngen',
                                                                 'seed','nfeatures','mode')]
  pre_refMatrix_cibersortx_args = arg_list[names(arg_list) %in% c('downsample','downsample_ratio')]

  benchmarking_obj = list()

  if (!is.null(bulkSimulator_methods)) {

    if('favilaco' %in% bulkSimulator_methods){
      stop('favilaco bulk simulation methods is no longer supported since deconvBenchmarking v0.1.1')
    }

    supported_methods <- c('homo', 'semi', 'heter', 'heter_sampleIDfree', 'immunedeconv', 'SCDC')
    if (any(!bulkSimulator_methods %in% supported_methods)) {
      stop(
        "Some bulkSimulation methods are not supported. Please provide only valid names. Supported methods are: ",
        paste(supported_methods, collapse = ", ")
      )
    }

    if (is.null(colnames_of_cellType) || is.na(colnames_of_cellType)) {
      stop('Please provide colnames_of_cellType when calling bulk simulation.')
    }

    if (!fixed_cell_type %in% unique(scMeta[,colnames_of_cellType])) {
      stop('Please make sure the provided "fixed_cell_type" is a cell type from scMeta.')
    }

    if (any(c('heter', 'semi') %in% bulkSimulator_methods)) {
      if (is.na(heter_cell_type)) {
        stop('Please provide "heter_cell_type" when using semi/heter simulation.')
      } else {
        if (!heter_cell_type %in% unique(scMeta[,colnames_of_cellType])) {
          stop('Please make sure the provided "heter_cell_type" is a cell type from scMeta.')
        }
      }
    }

    if('heter_sampleIDfree' %in% bulkSimulator_methods){

      if (!requireNamespace("scran", quietly = TRUE)){
        stop('Please make sure you have scran installed to run heter_sampleIDfree bulk simulation')
      }
    }

    if('immunedeconv' %in% bulkSimulator_methods){
      if (!requireNamespace("immunedeconv", quietly = TRUE)){
        stop('Please make sure you have immunedeconv installed to run immunedeconv bulk simulation')
      }
    }

    if('SCDC' %in% bulkSimulator_methods){
      if (!requireNamespace("SCDC", quietly = TRUE)){
        stop('Please make sure you have SCDC installed to run SCDC bulk simulation')
      }
    }

  }

  if('scran' %in% refMarkers_methods || 'scran' %in% refMatrix_methods){
    if (!requireNamespace("scran", quietly = TRUE)){
      stop('Please make sure you have scran installed to run refMarkers_scran')
    }

    if (!requireNamespace("BayesPrism", quietly = TRUE)){
      stop('Please make sure you have BayesPrism installed to run refMarkers_scran')
    }
  }

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
                                                       colnames_of_sample = colnames_of_sample,
                                                       simulated_frac = simulated_frac,
                                                       heter_cell_type  = heter_cell_type,
                                                       ncells_perSample=ncells_perSample,
                                                       min_chunkSize=min_chunkSize,
                                                       use_chunk = use_chunk,
                                                       colnames_of_subcluster = colnames_of_subcluster,
                                                       export_cellUsage = export_cellUsage,
                                                       nbulk = nbulk,
                                                       n.core = n.core),
                                                  bulkSimulator_args))

    cat('\n')

  }else{
    bulk_simulation_obj = list()
  }

  if(include_tcga == T){
    message('"include_tcga" option is no longer supported since deconvBenchmarking v0.1.1, please refer to build_tcga_obj() function to manually create a tcga object')
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

  if(!is.null(refMatrix_methods)){
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
#' @param marker_based_methods a character vector indicating which methods to use for marker_based deconvolution.
#'    Use \code{\link{list_deconv_marker}} to check for available method names and suggested packages associated them.
#'    Make sure you have the required packages installed to use these methods. Set to NULL if no bulk simulation is needed. Set to NULL if no marker_based deconvolution is needed.
#' @param regression_based_methods a character vector indicating which methods to use for regression_based deconvolution.
#'    Use \code{\link{list_deconv_regression}} to check for available method names and suggested packages associated them.
#'    Make sure you have the required packages installed to use these methods. Set to NULL if no bulk simulation is needed. Set to NULL if no regression_based deconvolution is needed.
#' @param cibersort_path path to CIBERSORT.R. This argument is required for 'cibersort' regression method
#' @param scExpr Single-cell expression matrix used to simulate bulk data the benchmarking_obj, with genes in rows and cells in columns. This argument is required for 'MuSiC'
#'    as they rely on a reference single-cell profile in their methods. Note that only training cells in the benchmarking_obj will be used to create a single-cell reference for 'MuSiC'
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr. This argument is required for 'MuSiC'
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta. This argument is required for 'MuSiC'
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. This argument is required for 'MuSiC'
#' @param refFree_methods a character vector indicating which methods to use for reference free deconvolution. Use \code{\link{list_deconv_refFree}} to check for available method names
#'    and suggested packages associated them. Make sure you have the required packages installed to use these methods.Set to NULL if no bulk simulation is needed..
#' @param k argument for reference free methods: number of cell types in bulk expression. If set to NA, will be default to true number of cell types in bulk expression
#' @param immunedeconv_methods a character vector indicating which methods to use from immunedeconv package. Available methods include:
#'    'xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'. Make sure you have the suggested package'immunedeconv' installed before using immunedeconv methods.
#'    Set to NULL if no immunedeconv_methods is needed.
#' @param tcga_abbreviation a character string indicating tcga-abbreviation for the bulk data to be deconvoluted, for example 'skcm'. Required for 'timer' and 'consensus_tme'.
#' @param n.core number of cores to use for parallel programming. Default = 1
#' @param ... additional arguments to be passed to the following functions: deconv_marker(), deconv_regression(), deconv_refFree(), and immunedeconv::deconvolute()
#'
#' @return a list of deconvolution results for each bulk expression in the benchmarking_obj
#' @export
#'
#' @examples
#' \dontrun{
#' # first create a benchmarking_obj using benchmarking_init():
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
#'                   # arguments for bulk simulation: select bulk simulation methods
#'                   bulkSimulator_methods = c('homo','semi','heter','heter_sampleIDfree'),
#'
#'                   # argument required for semi/heter_sampleIDfree bulk simulation method
#'                   heter_cell_type = 'malignant',
#'
#'                   # general simulation parameters
#'                   ncells_perSample = 500,
#'                   min_chunkSize = 20,
#'                   use_chunk = 'random',
#'
#'                   # parameters for adjusting heterogeneity in sample ID-free bulk simulation
#'                   dirichlet_cs_par = 0.1,
#'                   min.subcluster.size = 20,
#'                   max.num.cs = NA,
#'
#'                   # argument for marker constructions
#'                   refMarkers_methods = c('limma','scran'),
#'
#'                   # arguments for signature matrics construction
#'                   refMatrix_methods = c('raw','limma','scran'),
#'
#'
#'                   # export files for autogeneS and cibersortx
#'                   create_autogeneS_input = T,
#'                   create_cibersortx_input = T,
#'
#'                   n.core = 4
#'                   )
#'
#' # perform deconvolution on the benchmarking_obj using all categories of deconvolution methods, which includes:
#' # 1) marker-based methods: 'firstPC','gsva','debCAM','TOAST'
#' # 2) regression-based methods: 'nnls','cibersort','MuSiC','wRLM','RPC'
#' # 3) reference free methods: 'linseed','debCAM'
#' # 4) other deconvolution methods from immunedeconv package:
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
#'                     # arguments for other deconvolution methods
#'                     immunedeconv_methods = c('xcell','mcp_counter','epic','quantiseq','timer','abis','consensus_tme','estimate'),
#'                     tcga_abbreviation = 'SKCM', # arguments required for 'timer' and 'consensus_tme'
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
#'                     immunedeconv_methods = NULL,
#'
#'                     n.core = 4)
#'
#' }
benchmarking_deconv = function(benchmarking_obj,
                               # arguments for deconv_marker()
                               marker_based_methods = c('firstPC','gsva','debCAM','TOAST'),

                               # arguments for deconv_regression()
                               regression_based_methods =  c('nnls','cibersort','MuSiC','wRLM','RPC'),
                               cibersort_path = NULL,
                               scExpr = NULL, scMeta = NULL, colnames_of_cellType = NA, colnames_of_sample = NA,

                               # arguments for deconv_refFree()
                               refFree_methods = c('linseed','debCAM'),
                               k = NA,

                               # arguments for methods from immunedeconv package
                               immunedeconv_methods = c('xcell','mcp_counter','epic','quantiseq','timer','abis',
                                                 'consensus_tme','estimate'),
                               tcga_abbreviation = NA,
                               n.core = 1,
                               ...

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
    if(!require(immunedeconv, quietly = T)){
      message("The suggested package immunedeconv is not installed, will install it automatically")
      remotes::install_github("omnideconv/immunedeconv")
      require(immunedeconv)
    }else{
      require(immunedeconv, quietly = T)
    }

  }

  if('consensus_tme' %in% immunedeconv_methods | 'timer' %in% immunedeconv_methods){
    if(is.na(tcga_abbreviation)){
      stop('please provide tcga_abbreviation to run consensus_tme and/or timer')
    }
  }

  arg_list = list(...)
  deconv_marker_args = arg_list[names(arg_list) %in% c('alpha','sigma','epsilon','maxIter')]
  deconv_regression_args = arg_list[names(arg_list) %in% c('QN_cibersort','skip_raw_cibersort','normalize_MuSiC',
                                                           'skip_raw_wRLM','weight_wRLM','intercept_wRLM','scale_wRLM',
                                                           'QN_wRLM','skip_raw_RPC','maxit_RPC')]
  deconv_refFree_args = arg_list[names(arg_list) %in% c('corner.strategy','dim.rdc','thres.low','thres.high','cluster.method',
                                                        'cluster.num','MG.num.thres','lof.thres','quickhull','quick.select',
                                                        'sample.weight','appro3','generalNMF','cores')]

  immunedeconv_deconvolute_args = arg_list[names(arg_list) %in% c('tumor','arrays','column','rmgenes',
                                                                  'scale_mrna','expected_cell_types')]

  if(!is.null(scExpr) & !is.null(scMeta)){

    stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

    if(all(benchmarking_obj$splitting$training_cells %in% colnames(scExpr))){
      scExpr_train = scExpr[,benchmarking_obj$splitting$training_cells]
      scMeta_train = scMeta[benchmarking_obj$splitting$training_cells,]
    }else{
      scExpr_train = scExpr
      scMeta_train = scMeta
      warning('The provided scRNA profiles does not include all the training cells used in the splitting steps.
              The provided scRNA profile will be used for MuSiC')
    }
  }else{
    scExpr_train = NULL
    scMeta_train = NULL
  }

  for(bulk in names(benchmarking_obj$bulk_simulation_obj)){

    bulk_expr = benchmarking_obj[['bulk_simulation_obj']][[bulk]][['simulated_bulk']]
    bulk_deconvRes = list()

    if(!is.null(marker_based_methods)){
      marker_based_deconvRes = do.call(deconv_marker,c(list(methods = marker_based_methods,
                                                            bulk_expr = bulk_expr,
                                                            marker_list = benchmarking_obj$marker_list,
                                                            n.core = n.core),
                                                       deconv_marker_args))

      bulk_deconvRes = c(bulk_deconvRes,marker_based_deconvRes)
    }

    if(!is.null(regression_based_methods)){

      regression_based_deconvRes = do.call(deconv_regression,c(list(methods = regression_based_methods,
                                                                    bulk_expr = bulk_expr,
                                                                    sigMatrix_list = benchmarking_obj$sigMarix_list,
                                                                    cibersort_path = cibersort_path,
                                                                    scExpr = scExpr_train,
                                                                    scMeta = scMeta_train,
                                                                    colnames_of_cellType = colnames_of_cellType,
                                                                    colnames_of_sample = colnames_of_sample,
                                                                    n.core = n.core),
                                                               deconv_regression_args))

      bulk_deconvRes = c(bulk_deconvRes,regression_based_deconvRes)
    }

    if(!is.null(refFree_methods)){
      if(is.na(k)){
        k = ncol(benchmarking_obj[['bulk_simulation_obj']][[bulk]][['simulated_frac']])
      }

      refFree_deconvRes = do.call(deconv_refFree,c(list(methods = refFree_methods,
                                                        bulk_expr = bulk_expr,
                                                        k = k),
                                                   deconv_refFree_args))

      bulk_deconvRes = c(bulk_deconvRes,refFree_deconvRes)
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

        immunedeconv_deconvRes[[method]] = do.call(immunedeconv::deconvolute,c(list(gene_expression = bulk_expr,
                                                                                    method = method,
                                                                                    indications = indications),
                                                                               immunedeconv_deconvolute_args)) %>% column_to_rownames('cell_type') %>% as.matrix()

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
#' @param nonzero_threshold The threshold percentage of non-zero count in the simulated fraction. Only cell types with a non-zero count percentage equal to
#'    or higher than this threshold will be included for comparing their deconvolution performance.
#'
#' @return a list of performance evaluation statistics including per-cell type correlation and RMSE values, specified for each simulated bulk expression
#' @export
#'
#' @examples
#' \dontrun{
#' # a complete pipeline to benchmarking deconvolution methods using bulk data simulated under various strategies
#'
#' # first create a benchmarking_obj using benchmarking_init():
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
#'                   # arguments for bulk simulation: select bulk simulation methods
#'                   bulkSimulator_methods = c('homo','semi','heter','heter_sampleIDfree'),
#'
#'                   # argument required for semi/heter_sampleIDfree bulk simulation method
#'                   heter_cell_type = 'malignant',
#'
#'                   # general simulation parameters
#'                   ncells_perSample = 500,
#'                   min_chunkSize = 20,
#'                   use_chunk = 'random',
#'
#'                   # parameters for adjusting heterogeneity in sample ID-free bulk simulation
#'                   dirichlet_cs_par = 0.1,
#'                   min.subcluster.size = 20,
#'                   max.num.cs = NA,
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
#'
#'                   # export files for autogeneS and cibersortx
#'                   create_autogeneS_input = T,
#'                   create_cibersortx_input = T,
#'
#'                   n.core = 4
#'                   )
#'
#' # perform deconvolution on the benchmarking_obj using all categories of deconvolution methods, which includes:
#' # 1) marker-based methods: 'firstPC','gsva','debCAM','TOAST'
#' # 2) regression-based methods: 'nnls','cibersort','MuSiC','wRLM','RPC'
#' # 3) reference free methods: 'linseed','debCAM'
#' # 4) other deconvolution methods from immunedeconv package:
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
benchmarking_evalu = function(deconvResults_obj, nonzero_threshold = 0.05){

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

    # Only cell types with a non-zero count ratio equal to or higher than 'nonzero_ratios' will be included for comparing their deconvolution performance.
    nonzero_ratios = colSums(Y>0, na.rm = T)/nrow(Y)
    ct_to_evalu = names(nonzero_ratios)[nonzero_ratios >= nonzero_threshold]
    Y_ = Y[,ct_to_evalu,drop = F]

    bulk_performance = list()

    # detailed per-method evaluation statistics
    deconvEvalu = list()
    for(method in names(deconvResults_obj[[bulk]]$deconvRes)){

      E = deconvResults_obj[[bulk]]$deconvRes[[method]]
      E = E[,order(colnames(E)),drop = F]

      if(ncol(Y) == ncol(E)){
        if(all(colnames(Y) == colnames(E))){
          # modify E according to nonzero_threshold if it contains the exact set of cell types as the true_frac
          E = E[,ct_to_evalu,drop = F]
        }else{
          # find cell types that exhibit the highest correlation with true_frac in Y_
          maxCorName = c()
          for(ct in colnames(Y_)){
            id = apply(cor(E,Y_[,ct],use = 'pairwise.complete.obs'),2,which.max)
            maxCorName = c(maxCorName, colnames(E)[id])
          }

          E = E[,maxCorName,drop = F]
          colnames(E) = make.unique(colnames(E))
        }
      }else if(all(colnames(E) %in% colnames(Y))){
        # consider marker based methods when specific cell types are excluded because no enough number of markers pass the minimum_n threshold
        # report these cell types as NA in the performance export
        E = E[,colnames(E) %in% ct_to_evalu, drop = F]
        missing_ct = ct_to_evalu[!ct_to_evalu %in% colnames(E)]

        if(length(missing_ct)>0){
          empty_matrix = matrix(NA,ncol = length(missing_ct),nrow = nrow(E),dimnames = list(rownames(E),missing_ct))
          E = cbind(E,empty_matrix)
          E = E[,order(colnames(E)),drop = F]
        }

      }else{
        # consider other deconvolution methods
        maxCorName = c()
        for(ct in colnames(Y_)){
          id = apply(cor(E,Y_[,ct],use = 'pairwise.complete.obs'),2,which.max)
          maxCorName = c(maxCorName, colnames(E)[id])
        }

        E = E[,maxCorName,drop = F]
        colnames(E) = make.unique(colnames(E))
      }

      m1 = gather(Y_ %>% as.data.frame(),cell_type,true_frac)
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
    rownames(bulk_performance$maxCor) = colnames(Y_)

    bulk_performance$RMSE = do.call(cbind, lapply(deconvEvalu,get_RMSE))
    rownames(bulk_performance$RMSE) = colnames(Y_)

    bulk_performance$deconvEvalu = deconvEvalu


    deconvPerformance[[bulk]] = bulk_performance
    message(paste('>>>>>>>>>>>>>>>>>>>>>>>>>>>>> finish performance evaluation for',bulk,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))

  }

  return(deconvPerformance)
}


#' Summarize deconvolution performance in a long table
#'
#' @param deconvPerformance a deconvPerformance object returned from benchmarking_evalu() function
#'
#' @return a dataframe summarizing the deconvolution performance per deconvolution method and simulation method
#' @export
gather_performance = function(deconvPerformance){
  quick_gather = function(name, deconvPerformance){
    deconvEvalu = deconvPerformance[[name]][['deconvEvalu']]
    extract_summ = function(method){
      l = deconvEvalu[[method]][['summ']] %>% as.data.frame()  %>% rownames_to_column('cell_type') %>% mutate(method = method)
      return(l)
    }
    out = do.call(rbind,lapply(names(deconvEvalu),extract_summ)) %>% mutate(class = name)
  }
  out = do.call(rbind,lapply(names(deconvPerformance),quick_gather,deconvPerformance))
  return(out)
}

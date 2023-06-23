
create_train_test_splitting<-function(scMeta, colnames_of_cellType, training_ratio = 0.5){
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

  return(splitting_list)
}

#' Create an object for deconvolution benchmarking
#' @description This function takes scRNA profile as input and generate an object intended for future deconvolution benchmarking.
#'    The function performs the following steps: (1) It divides the cells into training and testing cells;
#'    (2) the training cells are utilized to generate reference profiles, such as markers and signature matrices;
#'    (3) the testing cells are used to generate simulated bulk expression, which is then employed for deconvolution purposes;
#'    (4) additionally, this function offers the flexibility to include a TCGA cohort as part of the object for future deconvolution benchmarking
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#'
#' @param fixed_cell_type argument for fracSimulator_Beta() function: a character denotes the target cell type for which we strive to faithfully preserve its distribution.
#'    It is recommended to set this parameter to the name of the malignant cell types. If left undefined, the function will automatically
#'    select the most abundant cell type as 'fixed_cell_type'.
#' @param min.frac fracSimulator_Beta() argument: minimum fraction in the simulated fraction, values below this threshold will be set to zero. Default = 0.01
#' @param showFractionPlot fracSimulator_Beta() argument: a logical variable determining whether to display simulated fraction distribution for the fixed_cell_type
#'
#' @param training_ratio ratio of training cells. Default = 0.5
#'
#' @param bulkSimulator_methods bulkSimulator() argument: a character vector indicating which bulk simulation methods to use. Use list_bulkSimulator() to check for available method names.
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_cellState column name that corresponds to cellState in scMeta, where cellState contains sub-clustering information for each cellType.  This is an argument required for 'heter_sampleIDfree' method only.
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. Required for 'semi', 'heter', and 'SCDC' methods.
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns. This argument if required for 'homo', 'semi', 'heter', 'heter_sampleIDfree', 'immunedeconv' methods
#' @param heter_cell_type name of the cell_type to maintain the highest level of heterogeneity. It is recommended to set this parameter to the name of the malignant cell types.
#'     This argument is required for 'semi' and 'heter_sampleIDfree' methods
#' @param ncells_perSample number of cells to aggregate for each simulated bulk sample. This is an argument required for 'homo', 'semi', 'favilaco', 'immunedeconv' and 'SCDC' methods
#' @param min_chunkSize minimum number of cells to aggregate to construct a given cell-type component in the simulated bulk. This is an argument required for 'semi' and 'heter' methods
#' @param use_chunk a character indicating which cells to pool together for the a given cell_type. Default='all' other options include 'random'.
#'    When use_chunk = 'all', use all the cells belonging to the same patient for a given cell type to generate the certain cell type component in the simulated bulk;
#'    when use_chunk = 'random', randomly select 50-100% of the cells belonging to the same patient for a given cell type. This is an argument required for 'semi' and 'heter' methods
#' @param dirichlet_cs_par a numeric value determine the dispersion level of the simulated fractions. With lower value indicating higher dispersion level. Default = 0.1. This is an argument required for 'heter_sampleIDfree' method.
#' @param min.percentage minimum percentage of cellType fraction to generate in fraction simulation. Default = 1. This argument is only required for 'favilaco'
#' @param max.percentage maximum percentage of cellType fraction to generate in fraction simulation. Default = 99. This argument is only required for 'favilaco'
#' @param nbulk number of simulated bulk samples. This argument is required for 'favilaco' and 'SCDC' methods.
#' @param seed a seed value for 'favilaco' method. Default = 24
#' @param use_simulated_frac_as_prop_mat a logical variable to determine whether to use 'simulated_frac' as 'prop_mat' as input for 'SCDC' method
#' @param disease indicate the health condition of subjects. This argument is only required for 'SCDC' method.
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used. This argument is only required for 'SCDC' method.
#' @param samplewithRep logical, randomly sample single cells with replacement. Default is T. This argument is only required for 'SCDC' method.
#' @param refMarkers_methods a character vector specifying the desired methods for generating cell-type specific markers. Use list_refMarkers() to check for available method names.
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes.
#' @param minimum_n minimum number of marker genes for a cell-type
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#' @param sigMatrixList a list of signature matices to derive markers from. This argument is required for 'sigMatrixList' method
#' @param refMatrix_methods a character vector specifying the desired methods for generating signature matrices. Use list_refMarix() to check for available method names.
#' @param markerList a list of pre-calculated cell-type marker list, where the first level of this list represents different methods used to derive cell-type specific markers, and the second level
#'    comprises the actual cell-type specific genes identified by each method. This argument is required for 'markerList' method
#' @param include_tcga a logicial variable determining whether to include tcga in the output object
#' @param tcga_abbreviation a character indicating tcga abbreviation for the tcga cohort to include, for example 'SKCM'
#' @param purity_method tumor purity estimation method that is utilized as a means of estimating the malignant proportion within the exported object for TCGA expression data.
#'    Available methods include 'ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC' and 'CPE', Default = 'CPE'
#' @param create_autogeneS_input a logical variable determine whether to create input data for autogeneS, which is a python based approach to construct signature matrix
#' @param max.spec_cutoff_for_autogeneS specificity score threshold to select for hv_genes. Default = 0.5
#' @param autogeneS_input_file_name desired file name to save the processed file
#' @param display_autogeneS_command a logical variable indicating whether to display the command lines to run autogeneS. By pasting the generated code into the command line, an external Python file will be executed,
#'    which will return the generated signature matrices in the 'autogeneS_output' folder.
#' @param ngen a numeric variable indicating number of generations used in autogeneS
#' @param seed_autogeneS autogeneS seed argument. Default = 0
#' @param nfeatures autogeneS nfeatures argument. Default = 400
#' @param mode autogeneS mode argument. Default = 'fixed'
#' @param create_cibersortx_input a logical variable determine whether to create input data for cibersortx, which is a web-server to construct signature matrix
#' @param downsample a logical variable indicating whether to downsample the processed file of not. This argument is useful when the processed file exceeds the storage limit for the cibersortx server
#' @param downsample_ratio a numeric value indicating downsampling ratio. This argument determines the fraction of the original size that will be retained in the downsampling process
#' @param cibersortx_input_file_name desired file name to save the processed file
#' @param n.core number of cores to use for parallel programming. Default = 1
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
#'                   fixed_cell_type = 'malignant',
#'                   bulkSimulator_methods = c('homo', 'semi','heter','favilaco','immunedeconv','SCDC'),
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # argument for semi bulk simulation method
#'                   heter_cell_type = 'malignant',
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
#'                   purity_method = 'CPE',
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
#'                   training_ratio = 1,
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
#'                   purity_method = 'CPE'
#'                   )
#' }
benchmarking_init = function(scExpr,scMeta,

                             # arguments for fracSimulator_Beta()
                             fixed_cell_type = NA,
                             min.frac = 0.01,
                             showFractionPlot = T,

                             # arguments for create_train_test_splitting()
                             training_ratio = 0.5,

                             # arguments for bulkSimulator()
                             bulkSimulator_methods = NULL,
                             colnames_of_cellType = NA,
                             colnames_of_cellState = NA,
                             colnames_of_sample = NA,
                             simulated_frac = NULL,
                             heter_cell_type = NA,
                             ncells_perSample = 500,
                             min_chunkSize = 5,
                             use_chunk = "all",
                             dirichlet_cs_par = 0.1,
                             min.percentage = 1,
                             max.percentage = 99,
                             nbulk = 100,
                             seed = 24,
                             use_simulated_frac_as_prop_mat = FALSE,
                             disease = NULL,
                             ct.sub = NULL,
                             samplewithRep = TRUE,

                             # arguments for refMarkers
                             refMarkers_methods = c('limma','scran'),
                             hv_genes = NULL,
                             log2FC = 2,
                             log2FC_flexible = 1,
                             minimum_n = 15,
                             maximum_n = 50,
                             max.spec_cutoff_for_DE = 0.3,
                             sigMatrixList = NULL,

                             # arguments for refMatrix
                             refMatrix_methods = c('raw','limma','scran'),
                             markerList = NULL,

                             # arguments for including TCGA as part of evaluation
                             include_tcga = F,
                             tcga_abbreviation = NA,
                             purity_method = 'CPE',

                             # arguments for pre_refMatrix_autogeneS
                             create_autogeneS_input = F,
                             max.spec_cutoff_for_autogeneS = 0.5,
                             autogeneS_input_file_name = NULL,
                             display_autogeneS_command = T,
                             ngen = 5000,
                             seed_autogeneS = 0,
                             nfeatures = 400,
                             mode = "fixed",

                             # arguments for pre_refMatrix_cibersortx
                             create_cibersortx_input = F,
                             downsample = F,
                             downsample_ratio = 0.2,
                             cibersortx_input_file_name = NULL,

                             n.core = 1
){
  benchmarking_obj = list()

  # split scRNA into training and testing
  splitting = create_train_test_splitting(scMeta, colnames_of_cellType, training_ratio)

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

  # generate a list of bulk simulation objects using testing scRNA
  if(!is.null(bulkSimulator_methods)){

    # Simulate realistic cell-type fractions from beta distribution
    if(is.null(simulated_frac)){
      simulated_frac = fracSimulator_Beta(scMeta,
                                          n = nbulk,
                                          colnames_of_sample,
                                          colnames_of_cellType,
                                          fixed_cell_type,
                                          min.frac,
                                          showFractionPlot)
    }

    bulk_simulation_obj = bulkSimulator(bulkSimulator_methods,
                                        scExpr_test,
                                        scMeta_test,
                                        colnames_of_cellType,
                                        colnames_of_cellState,
                                        colnames_of_sample,
                                        simulated_frac,
                                        heter_cell_type,
                                        ncells_perSample,
                                        min_chunkSize,
                                        use_chunk,
                                        dirichlet_cs_par,
                                        min.percentage,
                                        max.percentage,
                                        nbulk,
                                        seed,
                                        use_simulated_frac_as_prop_mat,
                                        disease,
                                        ct.sub,
                                        samplewithRep,
                                        n.core)
  }else{
    bulk_simulation_obj = list()
  }

  # include TCGA in the bulk object for further evaluation
  if(include_tcga == T){
    if(is.na(tcga_abbreviation)){
      stop('Please provide the TCGA cohort abbreviation, such as "SKCM", to specifically include this cohort in the bulk object')
    }
    tcga_abbreviation = stringr::str_to_upper(tcga_abbreviation)
    mapfile = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap'
    fpkmfile = paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-',tcga_abbreviation,'.htseq_fpkm.tsv.gz')
    to_use_genes = rownames(scExpr)

    message('downloading tcga expression from Xena browser')
    tcga_expr = xnea_fpkm2tpm(fpkmfile,mapfile,to_use_genes)
    purity_allMethods = tumor_purity(colnames(tcga_expr))

    purity = purity_allMethods[,purity_method][match(colnames(tcga_expr),rownames(purity_allMethods))]

    tcga_obj = list()
    tcga_obj$simulated_bulk = tcga_expr
    tcga_obj$simulated_frac = matrix(purity,ncol = 1,dimnames = list(NULL,purity_method))
    bulk_simulation_obj[[paste0('tcga_',tcga_abbreviation)]] = tcga_obj
  }

  # generate a list of cell-type specific markers from training scRNA
  if(!is.null(refMarkers_methods)){
    marker_list = refMarkers(refMarkers_methods,
                             scExpr_train,
                             cell_type_labels_train,
                             cell_state_labels_train,
                             hv_genes,
                             log2FC,
                             log2FC_flexible,
                             minimum_n,
                             maximum_n,
                             max.spec_cutoff_for_DE,
                             sigMatrixList)
  }else{
    marker_list = list()
  }

  # generate a list of signature matrices from training scRNA
  if(!is.null(refMarkers_methods)){
    if(all.equal(intersect(refMarkers_methods,refMatrix_methods),c('limma','scran'))){
      refMatrix_methods_remaining = refMatrix_methods[! refMatrix_methods %in% c('limma','scran')]
      sigMarix_list = refMatrix(methods = c(refMatrix_methods_remaining,'markerList'),
                                scExpr = scExpr_train,
                                cell_type_labels = cell_type_labels_train,
                                markerList = marker_list)
    }else{
      sigMarix_list = refMatrix(refMatrix_methods,
                                scExpr_train,
                                cell_type_labels_train,
                                cell_state_labels_train,
                                hv_genes,
                                log2FC,
                                log2FC_flexible,
                                minimum_n,
                                maximum_n,
                                max.spec_cutoff_for_DE,
                                markerList)
    }
  }else{
    sigMarix_list = list()
  }

  # prepare input for autogeneS using training scRNA
  if(create_autogeneS_input==T){
    pre_refMatrix_autogeneS(scExpr_train,
                            cell_type_labels_train,
                            hv_genes,
                            max.spec_cutoff_for_autogeneS,
                            autogeneS_input_file_name,
                            display_autogeneS_command,
                            ngen,
                            seed,
                            nfeatures,
                            mode)
  }

  # prepare input for cibersortx using training scRNA
  if(create_cibersortx_input ==T){
    pre_refMatrix_cibersortx(scExpr_train,
                             cell_type_labels_train,
                             downsample,
                             downsample_ratio,
                             scMeta,
                             colnames_of_sample,
                             colnames_of_cellType,
                             cibersortx_input_file_name)
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
#'                   fixed_cell_type = 'malignant',
#'                   bulkSimulator_methods = c('homo', 'semi','heter','favilaco','immunedeconv','SCDC'),
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # argument for semi bulk simulation method
#'                   heter_cell_type = 'malignant',
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
#'                   purity_method = 'CPE',
#'
#'                   n.core = 4)
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
#'                   fixed_cell_type = 'malignant',
#'                   bulkSimulator_methods = c('homo', 'semi','heter','favilaco','immunedeconv','SCDC'),
#'                   colnames_of_cellType = 'cell_type',
#'                   colnames_of_sample = 'sampleID',
#'
#'                   # argument for semi bulk simulation method
#'                   heter_cell_type = 'malignant',
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
#'                   purity_method = 'CPE',
#'
#'                   n.core = 4)
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




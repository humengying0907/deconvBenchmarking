
#' Construct signature matrix using average expression across cell types
#' @description This function calculates mean expression across cell types to form a signature matrix, without performing any gene filtering
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#'
#' @return a signature matrix with genes in rows and cell-types in columns
#' @export
refMatrix_raw = function(scExpr, cell_type_labels){
  C = build_ref_matrix(scExpr,cell_type_labels)
  return(C)
}

#' Construct signature matrix using markers derived from limma DE analysis
#' @description The function first calculates the mean expression across cells from the same cell types, forming a reference matrix C. The C matrix is then filtered
#'    based on markers genes identified from limma differential expression (DE) analysis.
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes.
#' @param minimum_n minimum number of marker genes for a cell-type
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#'
#' @return a signature matrix with genes in rows and cell-types in columns
#' @export
refMatrix_limma = function(scExpr,cell_type_labels, hv_genes=NULL,
                           log2FC=2, log2FC_flexible = 1, minimum_n=15,maximum_n=50,
                           max.spec_cutoff_for_DE = 0.3){
  limma_markers = refMarkers_limma(scExpr,cell_type_labels, hv_genes,
                                   log2FC, log2FC_flexible, minimum_n,maximum_n,
                                   max.spec_cutoff_for_DE)
  C = build_ref_matrix(scExpr,cell_type_labels)
  sig_genes = do.call(c,limma_markers)
  matrix_limma = C[sig_genes,]
  return(matrix_limma)
}


#' Construct signature matrix using markers derived from scran DE analysis
#' @description The function first calculates the mean expression across cells from the same cell types, forming a reference matrix C. The C matrix is then filtered
#'    based on markers genes identified from scran differential expression (DE) analysis.
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param cell_state_labels a character vector indicating cell-states of each cells in scExpr. Pairwise DE analysis is initially conducted between cell states and then merges for cell types.
#'    If set to NULL, it will default to the same values as the cell_type_labels.
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes.
#' @param minimum_n minimum number of marker genes for a cell-type
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#'
#' @return a signature matrix with genes in rows and cell-types in columns
#' @export
refMatrix_scran = function(scExpr,cell_type_labels,cell_state_labels = NULL,hv_genes = NULL,
                           log2FC=2, log2FC_flexible = 1, minimum_n = 15,maximum_n = 50,
                           max.spec_cutoff_for_DE = 0.3){
  scran_markers = refMarkers_scran(scExpr,cell_type_labels,cell_state_labels,hv_genes ,
                                   log2FC, log2FC_flexible, minimum_n,maximum_n,
                                   max.spec_cutoff_for_DE)
  C = build_ref_matrix(scExpr,cell_type_labels)
  sig_genes = do.call(c,scran_markers)
  matrix_scran = C[sig_genes,]
  return(matrix_scran)
}

#' Construct signature matrix with a pre-calculated marker list
#' @description  The function first calculates the mean expression across cells from the same cell types, forming a reference matrix C. The C matrix is then filtered
#'    based on pre-calculated markers
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param markerList a list of pre-calculated cell-type marker list, where the first level of this list represents different methods used
#'    to derive cell-type specific markers, and the second level comprises the actual cell-type specific genes identified by each method
#'
#' @return a list of signature matrix
#' @export
refMatrix_markerList = function(scExpr, cell_type_labels,markerList){
  C = build_ref_matrix(scExpr,cell_type_labels)
  out = list()
  for(name in names(markerList)){
    sig_genes = do.call(c,markerList[[name]])
    sig_genes = sig_genes[sig_genes %in% rownames(scExpr)]
    out = c(out,list(C[sig_genes,]))
  }
  names(out) = names(markerList)
  return(out)
}



#' Construct signature matrix from scRNA reference
#' @description Generate a list of signature matrices from scRNA reference with user-selected methods
#'
#' @param methods a character vector indicating which methods to use. Use list_refMarix() to check for available method names.
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param cell_state_labels a character vector indicating cell-states of each cells in scExpr. Pairwise DE analysis is initially conducted between cell states and then merges for cell types.
#'    If set to NULL, it will default to the same values as the cell_type_labels. This is an argument required for 'scran' method only
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes. This argument is required for method 'limma' and 'scran'
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold. This argument is required for method 'limma' and 'scran'
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes. This argument is required for method 'limma' and 'scran'
#' @param minimum_n minimum number of marker genes for a cell-type. This argument is required for method 'limma' and 'scran'
#' @param maximum_n maximum number of marker genes of a cell-type. This argument is required for method 'limma' and 'scran'
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3. This argument is required for method 'limma' and 'scran'
#' @param markerList a list of pre-calculated cell-type marker list, where the first level of this list represents different methods used to derive cell-type specific markers, and the second level
#'    comprises the actual cell-type specific genes identified by each method. This argument is required for 'markerList' method
#'
#' @return a list of signature matrices derived from scRNA reference
#' @export
#'
#' @examples
#' \dontrun{
#' refMatrix(methods = c('raw','markerList'), scExpr = scExpr, cell_type_labels = scMeta$cell_type, markerList = markers)
#' }
refMatrix = function(methods,
                     scExpr,cell_type_labels,
                     cell_state_labels = NULL,hv_genes = NULL,
                     log2FC = 2, log2FC_flexible = 1, minimum_n = 15,maximum_n = 50,
                     max.spec_cutoff_for_DE = 0.3, markerList = NULL){

  l = list_refMarix(show_description=T)
  sigMatrix_list = list()

  for(method in methods){
    description = l$description[l$method == method]
    message(description)

    switch(method,
           raw = {
             result = list('raw' = refMatrix_raw(scExpr, cell_type_labels))
           },
           limma = {
             result = list('limma' = refMatrix_limma(scExpr,cell_type_labels, hv_genes,
                                                      log2FC, log2FC_flexible, minimum_n,maximum_n,
                                                      max.spec_cutoff_for_DE))
           },
           scran = {
             result = list('scran' = refMatrix_scran(scExpr,cell_type_labels,cell_state_labels,hv_genes,
                                                      log2FC, log2FC_flexible, minimum_n,maximum_n,
                                                      max.spec_cutoff_for_DE))
           },
           markerList = {
             result = refMatrix_markerList(scExpr, cell_type_labels, markerList)
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_refMarix() to check for available methods"))
             result <- NULL
           })
    sigMatrix_list <- c(sigMatrix_list, result)
  }
  return(sigMatrix_list)
}


#' Prepare input for autogeneS
#' @description This function calculates the mean expression across cell types to form a signature matrix. It then filters the signature matrix by selecting highly variable genes
#'    to ensure computational efficiency; finally, it saves the resulting signature matrix in the 'autogeneS_input' folder for further signature matrix construction with autogeneS.
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec' threshold as hv_genes
#' @param max.spec_cutoff_for_autogeneS specificity score threshold to select for hv_genes. Default = 0.5
#' @param autogeneS_input_file_name desired file name to save the processed file
#' @param display_autogeneS_command a logical variable indicating whether to display the command lines to run autogeneS. By pasting the generated code into the command line, an external Python file will be executed,
#'    which will return the generated signature matrices in the 'autogeneS_output' folder.
#' @param ngen a numeric variable indicating number of generations used in autogeneS
#' @param seed autogeneS seed argument. Default = 0
#' @param nfeatures autogeneS nfeatures argument. Default = 400
#' @param mode autogeneS mode argument. Default = 'fixed'
#'
#' @return save the filtered signature matrix in the 'autogeneS_input' folder for further signature matrix construction with autogeneS
#' @export
pre_refMatrix_autogeneS <- function(scExpr, cell_type_labels, hv_genes = NULL, max.spec_cutoff_for_autogeneS = 0.5, autogeneS_input_file_name = NULL,
                                    display_autogeneS_command = T, ngen = 5000, seed = 0, nfeatures = 400, mode = 'fixed'){
  if(is.null(autogeneS_input_file_name)){
    autogeneS_input_file_name = 'centroids_sc_hv'
    warning('autogeneS_input_file_name is default to centroids_sc_hv.csv')
  }

  if(is.null(hv_genes)){
    max.spec =  compute.specificity(collapse(ref = scExpr %>% t(), labels = cell_type_labels))
    hv_genes = names(max.spec)[max.spec>max.spec_cutoff_for_autogeneS]
  }else{
    hv_genes = hv_genes[hv_genes %in% rownames(scExpr)]
  }

  cat(paste('select',length(hv_genes),'highly variable genes as input genes for autogeneS \n'))
  cat('\n')
  if(length(hv_genes) < 1000){
    warning('fewer than 1000 genes will be passed to autogeneS. Consider decreasing the value of max.spec to allow for the selection of a larger number of genes')
  }

  C = build_ref_matrix(scExpr,cell_type_labels)
  C = C[hv_genes,]

  if(display_autogeneS_command ==T){
    cat('paste the following code into the command line to construct signature matrix with autogeneS; refMatrix_autogeneS.py can be found at https://github.com/humengying0907/deconvBenchmarking/tree/main/scripts \n')
    print(paste('python refMatrix_autogeneS.py','--input_dir autogeneS_input','--ngen',ngen,'--seed', seed,'--nfeatures', nfeatures,'--mode', mode))
  }

  if(file.exists('autogeneS_input')){
    write.table(C,file = paste0('autogeneS_input/',autogeneS_input_file_name,'.csv'),sep = ',')
  }else{
    dir.create('autogeneS_input')
    write.table(C,file = paste0('autogeneS_input/',autogeneS_input_file_name,'.csv'),sep = ',')
  }
}


#' Prepare input for cibersortx server
#' @description This function pre-processes the single-cell RNA profile and saves the processed file in the 'cibersortx_input' folder. The saved file is in the correct format
#'    to be used as input for constructing a signature matrix with the cibersortx online portal
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param downsample a logical variable indicating whether to downsample the processed file of not. This argument is useful when the processed file exceeds the storage limit for the cibersortx server
#' @param downsample_ratio a numeric value indicating downsampling ratio. This argument determines the fraction of the original size that will be retained in the downsampling process
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. Required for 'semi', 'heter', and 'SCDC' methods.
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param cibersortx_input_file_name desired file name to save the processed file
#'
#' @return save the processed file in 'cibersortx_input' folder
#' @export
pre_refMatrix_cibersortx <- function(scExpr, cell_type_labels,
                                     downsample = F, downsample_ratio = 0.2,
                                     scMeta = NULL,
                                     colnames_of_sample = NA,
                                     colnames_of_cellType = NA,
                                     cibersortx_input_file_name = NULL){

  if(is.null(cibersortx_input_file_name)){
    cibersortx_input_file_name = 'cibersortx_input'
    warning('cibersortx_input_file_name is default to cibersortx_input.csv')
  }

  if(file.exists('cibersortx_input') ==F){
    dir.create('cibersortx_input')
  }

  if(downsample ==T){
    if(any(sapply(list(scMeta,colnames_of_sample,colnames_of_cellType), is.null))){
      warning('scMeta info is not available, complete random downsampling will be conducted')
      downsampled_id = sample(1:length(cell_type_labels),floor(length(cell_type_labels) * downsample_ratio))
    }else{
      downsampling<-function(scMeta,colnames_of_sample,colnames_of_cellType,downsample_ratio){
        scMeta_summary=scMeta %>% group_by_(colnames_of_sample, colnames_of_cellType) %>%
          summarise(n=n()) %>%
          mutate(freq=n/sum(n)) %>%group_by_(colnames_of_cellType) %>% summarise(sum=sum(n))
        expected_n=floor(nrow(scMeta)*downsample_ratio)
        perCellType_n=floor(expected_n/nrow(scMeta_summary))

        cellTypes=unique(scMeta[,colnames_of_cellType])
        downSampled_cells=c()

        for (cellType in cellTypes){
          scMeta_sub=scMeta[scMeta[,colnames_of_cellType]==cellType,]
          size = ifelse(nrow(scMeta_sub)>=perCellType_n,perCellType_n,nrow(scMeta_sub))
          # print(paste(cellType,size))
          downSample_cells_sub=sample(rownames(scMeta_sub),size = size)
          downSampled_cells=c(downSampled_cells,downSample_cells_sub)
        }

        remaining_n=expected_n-length(downSampled_cells)
        if (remaining_n>0){
          cells_candidate=rownames(scMeta)[!rownames(scMeta) %in% downSampled_cells]
          downSample_cells_sub=sample(cells_candidate,remaining_n)
          downSampled_cells=c(downSampled_cells,downSample_cells_sub)
        }
        downSampled_cells
      }

      downsampled_cell = downsampling(scMeta,colnames_of_sample,colnames_of_cellType,downsample_ratio)
      downsampled_id = match(downsampled_cell,colnames(scExpr))
    }
    sc_raw = scExpr[,downsampled_id]
    colnames(sc_raw)=cell_type_labels[downsampled_id]
    sc_raw=cbind(data.frame(GeneSymbol=rownames(sc_raw)),sc_raw)
    write.table(sc_raw,file = paste0('cibersortx_input/',cibersortx_input_file_name,'.txt'),sep = '\t',row.names = F)
  }else{
    sc_raw = scExpr
    colnames(scExpr) = cell_type_labels
    sc_raw=cbind(data.frame(GeneSymbol=rownames(sc_raw)),sc_raw)
    write.table(sc_raw,file = paste0('cibersortx_input/',cibersortx_input_file_name,'.txt'),sep = '\t',row.names = F)
  }
}


##################### limma markers ################
get_limma_statistics<-function(scExpr,cell_type_labels,hv_genes){
  stopifnot(ncol(scExpr)==length(cell_type_labels))
  stopifnot(max(scExpr)>100)
  annotation=factor(cell_type_labels)
  design <- model.matrix(~0+annotation)
  colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
  cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
  colnames(cont.matrix) <- colnames(design)
  diag(cont.matrix) <- (ncol(design)-1)/ncol(design)

  fit <- limma::lmFit(log2(scExpr[hv_genes,]+1), design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, trend=TRUE)

  return(list(fit2,cont.matrix))
}

get_limma_markers<-function(DE_statistics,log2FC=2, minimum_n=15,maximum_n=50){
  fit2=DE_statistics[[1]]
  cont.matrix=DE_statistics[[2]]

  markers=marker.fc(fit2,cont.matrix,log2.threshold = log2FC)
  markers=markers[markers$log2FC>0,]

  DE_list<-list()
  DE_cellTypes = unique(markers$CT)

  DE_cellTypes = DE_cellTypes[order(DE_cellTypes)]
  DE_cellTypes = DE_cellTypes[table(markers$CT)>=minimum_n]

  if(length(unique(markers$CT)) > length(DE_cellTypes)){
    warning(paste0('only ', length(DE_cellTypes),'/',length(unique(markers$CT)),' cell types containing DE genes based on current log2FC and minimum_n parameters'))
  }


  for (i in 1:length(DE_cellTypes)){

    DE_genes = rownames(markers)[markers$CT==DE_cellTypes[i]]
    if(length(DE_genes)>maximum_n){
      cat(paste('For',DE_cellTypes[i],':',length(DE_genes),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
      DE_list[[i]] = DE_genes[1:maximum_n]
    }else{
      DE_list[[i]] = DE_genes
    }

    names(DE_list)[i] = DE_cellTypes[i]

  }

  return(DE_list)
}

#' Obtain cell-type specific markers with limma DE analysis
#' @description Obtain cell-type specific markers by by conducting differential expression (DE) analysis using the limma package, which utilizes
#'    one-against-rest comparison for each cell types and filters for marker genes by selecting DE genes that meet a specified logFC threshold
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param minimum_n minimum number of marker genes for a cell-type. If a cell type has fewer than 'minimum_n' genes passing the log2FC threshold, it will be excluded from the marker list
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#'
#' @return a list of cell-type marker genes
#' @export
refMarkers_limma <-function(scExpr, cell_type_labels, hv_genes = NULL,
                            log2FC = 2, minimum_n = 15,maximum_n = 50,
                            max.spec_cutoff_for_DE = 0.3){
  if(is.null(hv_genes)){
    max.spec =  compute.specificity(collapse(ref = scExpr %>% as.matrix() %>% Matrix::t(), labels = cell_type_labels))
    hv_genes = names(max.spec)[max.spec>max.spec_cutoff_for_DE]
  }else{
    hv_genes = hv_genes[hv_genes %in% rownames(scExpr)]
  }
  require(limma, quietly = T) %>% suppressMessages()

  limma_statistics = get_limma_statistics(scExpr, cell_type_labels, hv_genes)
  limma_markers = get_limma_markers(limma_statistics, log2FC, minimum_n, maximum_n)

  return(limma_markers)
}

################### scran markers ####################
get_scran_statistics<-function(scExpr,cell_type_labels,cell_state_labels,hv_genes){

  if (!requireNamespace("scran", quietly = TRUE)){
    stop('Please make sure you have BayesPrism installed to run get_scran_statistics')
  }

  stopifnot(ncol(scExpr)==length(cell_type_labels))
  stopifnot(ncol(scExpr)==length(cell_state_labels))
  stopifnot(max(scExpr)>100)
  scExpr = as.matrix(scExpr)

  diff.exp.stat <- BayesPrism::get.exp.stat(sc.dat=scExpr[hv_genes,] %>% Matrix::t(),# filter genes to reduce memory use
                                            cell.type.labels=cell_type_labels,
                                            cell.state.labels=cell_state_labels,
                                            psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                            cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                            n.cores=1 #number of threads
  )

  return(diff.exp.stat)
}

get_scran_markers <- function(DE_statistics, log2FC = 2, minimum_n = 15,maximum_n=50, order_by = 'pval'){
  DE_list<-list()
  DE_cellTypes = names(DE_statistics)
  n = c()

  # count number of DE genes that pass the log2FC threshold in each cell type
  for(i in 1:length(DE_cellTypes)){
    n[i] = sum(DE_statistics[[i]]>log2FC)
  }
  names(n) = DE_cellTypes

  DE_cellTypes = names(n)[n>=minimum_n]

  for(ct in DE_cellTypes){
    if(order_by =='pval'){
      sub = DE_statistics[[ct]] %>% dplyr::arrange(desc(-pval.up.min))
      sub = sub[sub$min.lfc >= log2FC,]
      if(nrow(sub) > maximum_n){
        DE_list[[ct]] = rownames(sub)[1:maximum_n]
        cat(paste('For',ct,':',nrow(sub),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
      }else{
        DE_list[[ct]] = rownames(sub)
      }
    }else if(order_by == 'min.lfc'){
      sub = DE_statistics[[ct]] %>% dplyr::arrange(desc(min.lfc))
      sub = sub[sub$min.lfc >= log2FC,]
      if(nrow(sub) > maximum_n){
        DE_list[[ct]] = rownames(sub)[1:maximum_n]
        cat(paste('For',ct,':',nrow(sub),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
      }else{
        DE_list[[ct]] = rownames(sub)
      }
    }else{
      stop('Invalid "order_by" argument provided, please choose from either "pval" or "min.lfc"')
    }
  }
  return(DE_list)
}

#' Obtain cell-type specific markers with scran DE analysis
#' @description Obtain cell-type specific markers by by conducting differential expression (DE) analysis using the scran package, which utilizes
#'    pairwise comparison between each cell-types and filters for marker genes by selecting DE genes that meet a specified logFC threshold
#'
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param cell_state_labels a character vector indicating cell-states of each cells in scExpr. Pairwise DE analysis is initially conducted between cell states and then merges for cell types.
#'    If set to NULL, it will default to the same values as the cell_type_labels.
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param minimum_n minimum number of marker genes for a cell-type. If a cell type has fewer than 'minimum_n' genes passing the log2FC threshold, it will be excluded from the marker list
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param order_by a character specifying the criteria for ordering DE genes. Available options include 'pval' and 'min.lfc'. When more than 'maximum_n' genes pass the logFC threshold,
#'    markers will be selected from the top 'maximum_n' genes ordered by this criterion. Default = 'pval'.
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#'
#' @return a list of cell-type marker genes
#' @export
refMarkers_scran <- function(scExpr,cell_type_labels,cell_state_labels = NULL,hv_genes = NULL,
                             log2FC = 2, minimum_n = 15,maximum_n = 50, order_by = 'pval',
                             max.spec_cutoff_for_DE = 0.3){

  if (!requireNamespace("scran", quietly = TRUE)){
    stop('Please make sure you have scran installed to run refMarkers_scran')
  }

  if (!requireNamespace("BayesPrism", quietly = TRUE)){
    stop('Please make sure you have BayesPrism installed to run refMarkers_scran')
  }

  require(scran, quietly = T)

  if(is.null(cell_state_labels)){
    cell_state_labels = cell_type_labels
  }
  if(is.null(hv_genes)){
    max.spec =  compute.specificity(collapse(ref = scExpr %>% as.matrix() %>% Matrix::t(), labels = cell_type_labels))
    hv_genes = names(max.spec)[max.spec>max.spec_cutoff_for_DE]
  }else{
    hv_genes = hv_genes[hv_genes %in% rownames(scExpr)]
  }

  require(BayesPrism, quietly = T) %>% suppressMessages()

  scran_statistics = get_scran_statistics(scExpr,cell_type_labels,cell_state_labels,hv_genes)
  scran_markers = get_scran_markers(scran_statistics, log2FC, minimum_n, maximum_n, order_by)
  return(scran_markers)
}

############### Seurat markers ###########
#' Obtain cell-type specific markers using Seurat::FindAllMarkers()
#'
#' @param scExpr Single-cell expression data, with genes in rows and samples in columns
#' @param cell_type_labels a vector indicating cell-type level annotations
#' @param seurat_marker_method Seurat FindAllMarkers() test.use parameter. Default = 'wilcox'
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param minimum_n minimum number of marker genes for a cell-type. If a cell type has fewer than 'minimum_n' genes passing the log2FC threshold, it will be excluded from the marker list
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    Set to NULL if not available and the function will automatically find hv genes using Seurat::FindVariableFeatures() function
#' @param n.HVG Number of features to select as top variable features: Seurat::FindVariableFeatures() function nfeatures parameter. Default to 3000
#' @param ... additional parameters pass to Seurat::FindAllMarkers()
#'
#' @return a list of cell-type specific markers
#' @export
#'
refMarkers_Seurat = function(scExpr,cell_type_labels,
                             seurat_marker_method = 'wilcox',log2FC = 2,
                             minimum_n=15,maximum_n=50 ,hv_genes = NULL, n.HVG = 3000, ...){
  require(Seurat)
  seurat_obj = CreateSeuratObject(counts = scExpr)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- ScaleData(seurat_obj,features = rownames(seurat_obj))
  Idents(seurat_obj) = cell_type_labels

  if(is.null(hv_genes)){
    seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = n.HVG)
    marker <- FindAllMarkers(seurat_obj, features = VariableFeatures(seurat_obj), only.pos = T,
                             test.use = seurat_marker_method, logfc.threshold = log2FC, ...)
  }else{
    marker <- FindAllMarkers(seurat_obj, features = hv_genes, only.pos = T,
                             test.use = seurat_marker_method, logfc.threshold = log2FC,...)
  }

  gc()

  marker = marker[order(marker$cluster,marker$avg_log2FC,decreasing = T),]

  DE_list<-list()
  DE_cellTypes = unique(marker$cluster)
  DE_cellTypes = DE_cellTypes[order(DE_cellTypes)]
  DE_cellTypes = DE_cellTypes[table(marker$cluster)>=minimum_n]

  if(length(unique(marker$cluster)) > length(DE_cellTypes)){
    warning(paste0('only ', length(DE_cellTypes),'/',length(unique(marker$cluster)),' cell types containing DE genes based on current log2FC and minimum_n parameters'))
  }


  for (i in 1:length(DE_cellTypes)){

    DE_genes = marker$gene[marker$cluster == DE_cellTypes[i]]

    if(length(DE_genes)>maximum_n){
      cat(paste('For',DE_cellTypes[i],':',length(DE_genes),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
      DE_list[[i]] = DE_genes[1:maximum_n]
    }else{
      DE_list[[i]] = DE_genes
    }

    names(DE_list)[i] = DE_cellTypes[i] %>% as.vector()

  }
  return(DE_list)
}

############## markers from a signature matrix ##########
#' Obtain cell-type specific markers from a list of signature matrices
#'
#' @param sigMatrixList a list of signature matrices
#' @param minimum_n minimum number of marker genes for a cell-type. If a cell type has fewer than 'minimum_n' marker genes, it will be excluded from the marker list
#' @param maximum_n maximum number of marker genes of a cell-type
#'
#' @return a list of markers obtained from different input signature matrices. Each list element is named after the corresponding signature matrix used to generate the markers.
#' @export
refMarkers_sigMatrixList = function(sigMatrixList,minimum_n = 15, maximum_n = 50){
  require(debCAM, quietly = T) %>% suppressMessages()

  marker_list = list()
  marker_list_names = c()
  for(name in names(sigMatrixList)){
    ref = sigMatrixList[[name]]
    fc = debCAM::MGstatistic(ref,colnames(ref))
    fc$gene = rownames(ref)

    marker_counts = table(fc$idx)
    ct_to_include = names(marker_counts)[marker_counts>=minimum_n]

    if(length(ct_to_include)<1){
      warning(paste('In signature matrix',name,', no cell type has markers exceeding the "minimum_n" threshold. Consider lower minimum_n.' ))
      next
    }

    extract_markers = function(cell_type){
      fc_sub = fc[fc$idx ==cell_type,]
      fc_sub = fc_sub[order(fc_sub$OVE.FC,decreasing = T),]

      if(nrow(fc_sub)>=maximum_n){
        m = fc_sub$gene[1:maximum_n]
      }else if(nrow(fc_sub) < maximum_n){
        m = fc_sub$gene
      }
      return(list(m))
    }

    v = do.call(c,lapply(ct_to_include,extract_markers))
    names(v) = ct_to_include

    marker_list = c(marker_list,list(v))
    marker_list_names = c(marker_list_names,name)
  }
  names(marker_list) = marker_list_names

  return(marker_list)
}

############### Seurat markers ###########
#' Obtain cell-type specific markers using Seurat::FindAllMarkers()
#'
#' @param scExpr Single-cell expression data, with genes in rows and samples in columns
#' @param cell_type_labels a vector indicating cell-type level annotations
#' @param seurat_marker_method Seurat FindAllMarkers() test.use parameter. Default = 'wilcox'
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param minimum_n minimum number of marker genes for a cell-type. If a cell type has fewer than 'minimum_n' genes passing the log2FC threshold, it will be excluded from the marker list
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    Set to NULL if not available and the function will automatically find hv genes using Seurat::FindVariableFeatures() function
#' @param n.HVG Number of features to select as top variable features: Seurat::FindVariableFeatures() function nfeatures parameter. Default to 3000
#' @param ... additional parameters pass to Seurat::FindAllMarkers()
#'
#' @return a list of cell-type specific markers
#' @export
#'
refMarkers_Seurat = function(scExpr,cell_type_labels,
                             seurat_marker_method = 'wilcox',log2FC = 2,
                             minimum_n=15,maximum_n=50 ,hv_genes = NULL, n.HVG = 3000, ...){
  seurat_obj = CreateSeuratObject(counts = scExpr)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- ScaleData(seurat_obj,features = rownames(seurat_obj))
  Idents(seurat_obj) = cell_type_labels

  if(is.null(hv_genes)){
    seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = n.HVG)
    marker <- FindAllMarkers(seurat_obj, features = VariableFeatures(seurat_obj), only.pos = T,
                             test.use = seurat_marker_method, logfc.threshold = log2FC, ...)
  }else{
    marker <- FindAllMarkers(seurat_obj, features = hv_genes, only.pos = T,
                             test.use = seurat_marker_method, logfc.threshold = log2FC,...)
  }

  gc()

  marker = marker[order(marker$cluster,marker$avg_log2FC,decreasing = T),]

  DE_list<-list()
  DE_cellTypes = unique(marker$cluster)
  DE_cellTypes = DE_cellTypes[order(DE_cellTypes)]
  DE_cellTypes = DE_cellTypes[table(marker$cluster)>=minimum_n]

  if(length(unique(marker$cluster)) > length(DE_cellTypes)){
    warning(paste0('only ', length(DE_cellTypes),'/',length(unique(marker$cluster)),' cell types containing DE genes based on current log2FC and minimum_n parameters'))
  }


  for (i in 1:length(DE_cellTypes)){

    DE_genes = marker$gene[marker$cluster == DE_cellTypes[i]]

    if(length(DE_genes)>maximum_n){
      cat(paste('For',DE_cellTypes[i],':',length(DE_genes),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
      DE_list[[i]] = DE_genes[1:maximum_n]
    }else{
      DE_list[[i]] = DE_genes
    }

    names(DE_list)[i] = DE_cellTypes[i] %>% as.vector()

  }
  return(DE_list)
}

############### get markers  ############
#' Obtain cell-type specific markers from scRNA reference
#' @description Generate a list of markers from scRNA reference with user-selected methods
#'
#' @param methods a character vector indicating which methods to use. Use list_refMarkers() to check for available method names.
#' @param scExpr Single-cell expression data to perform differential expression (DE) analysis, with genes in rows and samples in columns
#' @param cell_type_labels a character vector indicating cell-types of each cells in scExpr
#' @param cell_state_labels a character vector indicating cell-states of each cells in scExpr. Pairwise DE analysis is initially conducted between cell states and then merges for cell types.
#'    If set to NULL, it will default to the same values as the cell_type_labels. Required for 'scran' method
#' @param hv_genes a character vector containing the names of high-variable genes. scExpr will be pre-filtered based on the provided hv_genes to reduce computation time during the differential expression (DE) analysis.
#'    If set to NULL, the function will automatically select genes with specificity score passing 'max.spec_cutoff_for_DE' threshold as hv_genes
#' @param log2FC log fold change threshold to select marker genes. Marker genes will be limited to a maximum of 'maximum_n' genes among those that pass the 'log2FC' threshold.
#' @param minimum_n minimum number of marker genes for a cell-type. If a cell type has fewer than 'minimum_n' marker genes, it will be excluded from the marker list.
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param order_by a character specifying the criteria for ordering DE genes from scran DE analysis. Available options include 'pval' and 'min.lfc'. When more than 'maximum_n' genes pass the logFC threshold,
#'    markers will be selected from the top 'maximum_n' genes ordered by this criterion. Default = 'pval'.
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#' @param sigMatrixList a list of signature matrices to derive markers from. This argument is required for 'sigMatrixList' method
#'
#' @return a list of cell-type specific markers generated from use-defined methods
#' @export
#'
#' @examples
#' \dontrun{
#' # Obtain cell-type specific markers with DE analysis
#' refMarkers(methods = c('limma','scran'), scExpr = scExpr, cell_type_labels = scMeta$cell_type)
#'
#' # Obtain cell-type specific markers from a list of signature matrices
#' # In real practice, these signature matrices can be generated using various methods outside of the provided DE-based methods, such as from the cibersortx server or using autogeneS.
#' # For more details on how to prepare the scRNA data for cibersortx and autogeneS, use the following command:
#' help("pre_refMatrix_autogeneS")
#' help("pre_refMatrix_cibersortx")
#' refMarkers(method = 'sigMatrixList', sigMatrixList = sigMatrixList, maximum_n = 50)
#' }
refMarkers = function(methods,
                      scExpr = NULL,cell_type_labels = NULL,cell_state_labels = NULL,hv_genes = NULL,
                      log2FC=2, minimum_n=15, maximum_n=50, order_by = 'pval',
                      max.spec_cutoff_for_DE = 0.3,
                      sigMatrixList = NULL){
  if('scran' %in% methods){
    if (!requireNamespace("scran", quietly = TRUE)){
      stop('Please make sure you have scran installed to run refMarkers_scran')
    }

    if (!requireNamespace("BayesPrism", quietly = TRUE)){
      stop('Please make sure you have BayesPrism installed to run refMarkers_scran')
    }
  }

  if('sigMatrixList' %in% methods){
    if(is.null(sigMatrixList)){
      stop('please provide the sigMatrixList argument to extract markers from')
    }
  }

  if(any(methods %in% c('limma','scran'))){
    if(is.null(hv_genes)){
      max.spec =  compute.specificity(collapse(ref = scExpr %>% as.matrix() %>% Matrix::t(), labels = cell_type_labels))
      hv_genes = names(max.spec)[max.spec>max.spec_cutoff_for_DE]
    }else{
      hv_genes = hv_genes[hv_genes %in% rownames(scExpr)]
    }
  }

  l = list_refMarkers(show_description=T)

  m = list()

  for(method in methods){
    description = l$description[l$method == method]
    message(description)

    switch(method,
           limma = {
             result = refMarkers_limma(scExpr, cell_type_labels, hv_genes,
                                       log2FC, minimum_n, maximum_n,
                                       max.spec_cutoff_for_DE)
             result = list('limma' = result)
           },
           scran = {
             result = refMarkers_scran(scExpr, cell_type_labels, cell_state_labels, hv_genes,
                                       log2FC, minimum_n, maximum_n, order_by,
                                       max.spec_cutoff_for_DE)
             result = list('scran' = result)
           },
           sigMatrixList = {
             result = refMarkers_sigMatrixList(sigMatrixList,minimum_n,maximum_n)
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_refMarkers() to check for available methods"))
             result <- NULL
           })
    m <- c(m, result)
  }
  return(m)
}





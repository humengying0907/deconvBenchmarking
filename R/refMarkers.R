
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

get_limma_markers<-function(DE_statistics,log2FC=2, log2FC_flexible = 1, minimum_n=15,maximum_n=50){
  fit2=DE_statistics[[1]]
  cont.matrix=DE_statistics[[2]]

  markers_flexible=marker.fc(fit2,cont.matrix,log2.threshold = log2FC_flexible)
  markers_flexible=markers_flexible[markers_flexible$log2FC>0,]

  markers=marker.fc(fit2,cont.matrix,log2.threshold = log2FC,)
  markers=markers[markers$log2FC>0,]

  DE_list<-list()
  DE_cellTypes=unique(markers$CT)

  for (i in 1:length(DE_cellTypes)){
    first_try=rownames(markers)[markers$CT==DE_cellTypes[i]]
    second_try=rownames(markers_flexible)[markers_flexible$CT==DE_cellTypes[i]]

    if(length(first_try)>=minimum_n){
      if(length(first_try)>maximum_n){
        cat(paste('For',DE_cellTypes[i],':',length(first_try),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
        DE_list[[i]]=first_try[1:maximum_n]
      }else{
        DE_list[[i]]=first_try
      }
    }else if (length(second_try)>=minimum_n){
      DE_list[[i]]=second_try[1:minimum_n]
      cat(paste('For',DE_cellTypes[i],': not enough genes passing log2FC threshod, choose first',minimum_n,'genes that passed log2FC =',log2FC_flexible,'threshold \n'))
    }else{
      DE_list[[i]]=second_try
      warning(paste('For',DE_cellTypes[i],': not enough genes passing log2FC =',log2FC_flexible, 'threshold, choose all genes passing log2FC =',log2FC_flexible,'threshold \n'))
    }

    names(DE_list)[i]=DE_cellTypes[i]
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
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes.
#' @param minimum_n minimum number of marker genes for a cell-type
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#'
#' @return a list of cell-type marker genes
#' @export
refMarkers_limma <-function(scExpr, cell_type_labels, hv_genes=NULL,
                            log2FC=2, log2FC_flexible = 1, minimum_n=15,maximum_n=50,
                            max.spec_cutoff_for_DE = 0.3){
  if(is.null(hv_genes)){
    max.spec =  compute.specificity(collapse(ref = scExpr %>% t(), labels = cell_type_labels))
    hv_genes = names(max.spec)[max.spec>max.spec_cutoff_for_DE]
  }else{
    hv_genes = hv_genes[hv_genes %in% rownames(scExpr)]
  }
  require(limma)
  limma_statistics = get_limma_statistics(scExpr,cell_type_labels,hv_genes)
  limma_markers = get_limma_markers(limma_statistics, log2FC, log2FC_flexible, minimum_n,maximum_n)
  return(limma_markers)
}

################### scran markers ####################
get_scran_statistics<-function(scExpr,cell_type_labels,cell_state_labels,hv_genes){
  stopifnot(ncol(scExpr)==length(cell_type_labels))
  stopifnot(ncol(scExpr)==length(cell_state_labels))
  stopifnot(max(scExpr)>100)

  diff.exp.stat <- BayesPrism::get.exp.stat(sc.dat=scExpr[hv_genes,] %>% t(),# filter genes to reduce memory use
                                            cell.type.labels=cell_type_labels,
                                            cell.state.labels=cell_state_labels,
                                            psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                                            cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                                            n.cores=1 #number of threads
  )

  return(diff.exp.stat)
}

get_scran_markers <- function(DE_statistics,log2FC=2,log2FC_flexible =1, minimum_n=15,maximum_n=50){
  DE_list<-list()
  DE_cellTypes=names(DE_statistics)

  for (i in 1:length(DE_cellTypes)){
    sub=DE_statistics[[i]] %>% dplyr::arrange(desc(-pval.up.min))
    first_try=rownames(sub)[sub$min.lfc>log2FC]
    second_try=rownames(sub)[sub$min.lfc>log2FC_flexible]

    if(length(first_try)>=minimum_n){
      if(length(first_try)>maximum_n){
        cat(paste('For',DE_cellTypes[i],':',length(first_try),'genes passed log2FC threshold, will pick top',maximum_n, '\n'))
        DE_list[[i]]=first_try[1:maximum_n]
      }else{
        DE_list[[i]]=first_try
      }
    }
    else if (length(second_try)>=minimum_n){
      DE_list[[i]]=second_try[1:minimum_n]
      cat(paste('For',DE_cellTypes[i],': not enough genes passing log2FC threshod, choose first',minimum_n,'genes that passed log2FC =',log2FC_flexible,'threshold \n'))
    }else{
      DE_list[[i]]=second_try
      warning(paste('For',DE_cellTypes[i],': not enough genes passing log2FC =',log2FC_flexible, 'threshold, choose all genes passing log2FC =',log2FC_flexible,'threshold \n'))
    }

    names(DE_list)[i]=DE_cellTypes[i]

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
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes.
#' @param minimum_n minimum number of marker genes for a cell-type
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#'
#' @return a list of cell-type marker genes
#' @export
refMarkers_scran <- function(scExpr,cell_type_labels,cell_state_labels = NULL,hv_genes = NULL,
                             log2FC=2, log2FC_flexible =1, minimum_n=15,maximum_n=50,
                             max.spec_cutoff_for_DE = 0.3){
  if(is.null(cell_state_labels)){
    cell_state_labels = cell_type_labels
  }
  if(is.null(hv_genes)){
    max.spec =  compute.specificity(collapse(ref = scExpr %>% t(), labels = cell_type_labels))
    hv_genes = names(max.spec)[max.spec>max.spec_cutoff_for_DE]
  }else{
    hv_genes = hv_genes[hv_genes %in% rownames(scExpr)]
  }

  require(BayesPrism)
  require(scran)

  scran_statistics = get_scran_statistics(scExpr,cell_type_labels,cell_state_labels,hv_genes)
  scran_markers = get_scran_markers(scran_statistics, log2FC, log2FC_flexible, minimum_n, maximum_n)
  return(scran_markers)
}

############## markers from a signature matrix ##########
#' Obtain cell-type specific markers from a list of signature matrices
#'
#' @param sigMatrixList a list of signature matrices
#' @param maximum_n maximum number of marker genes of a cell-type
#'
#' @return a list of markers obtained from different input signature matrices. Each list element is named after the corresponding signature matrix used to generate the markers.
#' @export
refMarkers_sigMatrixList = function(sigMatrixList,maximum_n = 50){
  require(debCAM)

  marker_list = list()
  for(name in names(sigMatrixList)){
    ref = sigMatrixList[[name]]
    fc = debCAM::MGstatistic(ref,colnames(ref))
    fc$gene = rownames(ref)

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

    v = do.call(c,lapply(colnames(ref),extract_markers))
    names(v) = colnames(ref)

    marker_list = c(marker_list,list(v))
  }
  names(marker_list) = names(sigMatrixList)
  return(marker_list)
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
#' @param log2FC_flexible a flexible log fold change threshold to select marker genes. If there are fewer than 'minimum_n' genes that pass the 'log2FC_flexible' threshold,
#'    all the genes that pass the threshold will be considered as marker genes.
#' @param minimum_n minimum number of marker genes for a cell-type
#' @param maximum_n maximum number of marker genes of a cell-type
#' @param max.spec_cutoff_for_DE specificity score threshold to select for hv_genes. Default = 0.3
#' @param sigMatrixList a list of signature matices to derive markers from. This argument is required for 'sigMatrixList' method
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
                      log2FC=2, log2FC_flexible =1, minimum_n=15,maximum_n=50,
                      max.spec_cutoff_for_DE = 0.3,
                      sigMatrixList = NULL){
  if('sigMatrixList' %in% methods){
    if(is.null(sigMatrixList)){
      stop('please provide the sigMatrixList argument to extract markers from')
    }
  }

  if(any(methods %in% c('limma','scran'))){
    if(is.null(hv_genes)){
      max.spec =  compute.specificity(collapse(ref = scExpr %>% t(), labels = cell_type_labels))
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
                                       log2FC, log2FC_flexible, minimum_n, maximum_n,
                                       max.spec_cutoff_for_DE)
             result = list('limma' = result)
           },
           scran = {
             result = refMarkers_scran(scExpr,cell_type_labels,cell_state_labels,hv_genes,
                                       log2FC, log2FC_flexible, minimum_n,maximum_n,
                                       max.spec_cutoff_for_DE)
             result = list('scran' = result)
           },
           sigMatrixList = {
             result = refMarkers_sigMatrixList(sigMatrixList,maximum_n)
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_refMarkers() to check for available methods"))
             result <- NULL
           })
    m <- c(m, result)
  }
  return(m)
}





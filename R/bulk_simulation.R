
#################### functions for homo simulation ######################
setCells <- function(cf, setList, ncells_perSample=500) {
  x=cf*ncells_perSample
  sc <- c()
  for (cType in names(x)) {
    n <- ceiling(x[cType]) # ceiling ensure that  for cell types with very small fraction, we can still get at least one cell of this cell type
    if (n > 0) {
      repl <- ifelse(n > length(setList[[cType]]),TRUE,FALSE)
      sc <- c(sc,sample(setList[[cType]],size = n,replace = repl))  # randomly draw cells from one cell type
    }
  }
  return(sc)
}

setBulks<-function(ID,simulated_frac,CellNameList,scExpr,cell_type_labels){
  named_labels=cell_type_labels
  names(named_labels)=colnames(scExpr)

  simulated_frac_row=simulated_frac[ID,]
  p=CellNameList[[ID]]
  X=scExpr[,p]

  # to account for the imbalanced cells distribution during the ceiling step in setCells(), re-weigh the cell-type specific expression
  cell_type_labels_sub=named_labels[p] %>% unname()
  C=build_ref_matrix(X,cell_type_labels_sub)
  E=C %*% matrix(simulated_frac_row[colnames(C)],ncol = 1)
  E=apply(E,2,function(x)((x/sum(x))*1e+06))

  return(E)
}


#' Homogeneous bulk simulation
#' @description Generate bulk samples by aggregating single cells with pre-defined fractions
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns
#' @param ncells_perSample number of cells to pool together to generate a simulated bulk sample. Default 500
#' @param export_cellUsage a logical variable determining whether to export cell names used to generate the simulated bulk. Default = F
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export

bulkSimulator_homo<-function(scExpr,scMeta,
                             colnames_of_cellType = NA,
                             simulated_frac=NULL,ncells_perSample=500,
                             export_cellUsage = F,
                             n.core = 1){
  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.null(simulated_frac)){
    stop('please provide cell fraction matrix with pre-defined proportions for bulk simulation')
  }


  cell_type_labels=scMeta[,colnames_of_cellType]

  stopifnot(all(colnames(simulated_frac) %in% unique(cell_type_labels)))

  cellNames=colnames(scExpr)
  cellSet=list()
  for (ts in levels(factor(cell_type_labels))) {
    cellSet[[ts]] <- cellNames[cell_type_labels == ts]
  }

  n = nrow(simulated_frac)
  CellNameList <- lapply(apply(simulated_frac,1,setCells,cellSet,ncells_perSample),function(x) x)

  pboptions(type = "txt", style = 3, char = "=")
  D=do.call(cbind,pblapply(seq(1,n),setBulks,simulated_frac,CellNameList,scExpr,cell_type_labels,cl = n.core))
  colnames(D)=paste0('simulated',1:ncol(D))

  if(export_cellUsage){
    cellUsed = unlist(CellNameList)
    cell_usage = data.frame(cell = cellUsed,
                            bulk_id = paste0('simulated',rep(seq_along(CellNameList), lengths(CellNameList))))

    names(cell_type_labels) = rownames(scMeta)
    cell_usage$cell_type = cell_type_labels[match(cell_usage$cell,names(cell_type_labels))]

    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac,
                cell_usage = cell_usage))

  }else{
    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac))
  }
}

######################### functions for semiheter simulation #####################
setCountSemi <- function(ID,
                         simulated_frac,
                         scExpr,scMeta,
                         colnames_of_cellType,
                         colnames_of_sample,
                         heter_cell_type,
                         ncells_perSample=500,
                         min_chunkSize=10,
                         use_chunk = 'all',
                         export_cellUsage = F){

  simulated_frac_row=simulated_frac[ID,]

  scMeta_renamed = data.frame(row.names = rownames(scMeta),
                              sampleID = scMeta[,colnames_of_sample],
                              cell_type = scMeta[,colnames_of_cellType])

  fixed_metaSub=scMeta_renamed[scMeta_renamed$cell_type == heter_cell_type,]
  sampleIDs=unique(fixed_metaSub$sampleID)
  fixed_sample=sample(sampleIDs,1)
  fixed_cells=rownames(fixed_metaSub)[fixed_metaSub$sampleID==fixed_sample]

  # ensure that there's enough number of cells to construct the fix-cell-type specific expression
  scMeta_unselected = fixed_metaSub[fixed_metaSub$sampleID!= fixed_sample,]
  scMeta_unselected = meta_shuffle(scMeta_unselected)

  if(length(fixed_cells) < min_chunkSize){
    patch_n = min_chunkSize - length(fixed_cells)

    if(nrow(scMeta_unselected) < patch_n){
      fixed_cells = c(fixed_cells,rownames(scMeta_unselected))
    }else{
      fixed_cells = c(fixed_cells,rownames(scMeta_unselected)[1:patch_n])
    }
  }


  if(use_chunk=='all'){
    # use all cells
    fixed_cells = fixed_cells
  }else if(use_chunk =='random'){
    # randomly select 50%-100% cells
    percent=sample(seq(50,100),size = 1)
    fixed_cells=sample(fixed_cells,size = ceiling(length(fixed_cells)*percent*0.01))
  }


  # non-malignant cells drawn from all patients
  cellNames=rownames(scMeta)
  cell_type_labels=scMeta[,colnames_of_cellType]
  cellSet=list()
  for (ts in levels(factor(cell_type_labels))) {
    cellSet[[ts]] <- cellNames[cell_type_labels == ts]
  }

  cellSet_remaining=cellSet[names(cellSet)!=heter_cell_type]
  simulated_frac_row_remaining=simulated_frac_row[names(simulated_frac_row)!=heter_cell_type]

  x=simulated_frac_row_remaining*ncells_perSample
  sc <- c()
  for (cType in names(x)) {
    n <- ceiling(x[cType]) %>% as.numeric() # ceiling ensure that even for cell types with low fraction, there's still at least one cell
    if (n > 0) {
      repl <- ifelse(n > length(cellSet_remaining[[cType]]),TRUE,FALSE)
      sc <- c(sc,sample(cellSet_remaining[[cType]],size = n,replace = repl))  # randomly draw cells from one cell type
    }
  }

  # final cells to use
  sc=c(sc,fixed_cells)

  # generate simulated expression by matrix multiplication
  cellType <- scMeta[sc,colnames_of_cellType]
  group = list()
  for(i in unique(cellType)){
    group[[i]] <- which(cellType %in% i)
  }

  X = scExpr[,sc]
  C = lapply(group,function(x) Matrix::rowMeans(X[,x,drop=F]))
  C = do.call(cbind, C)

  E=C %*% matrix(simulated_frac_row[colnames(C)],ncol = 1)
  E=apply(E,2,function(x)((x/sum(x))*1e+06))

  if(export_cellUsage){
    return(list(E = E,
                cellUsed = sc))
  }else{
    return(E)
  }
}


#' Semi-heterogeneous bulk simulation
#' @description Generate bulk samples by aggregating single cells, while maintaining heterogeneity for a certain cell type.
#'    Specifically, heter_cell_type is constrained to cells originating from a single patient. This requires that sampleID is pre-defined
#'    in the scRNA dataset
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param heter_cell_type cell_type to maintain heterogeneity across samples
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns
#' @param ncells_perSample number of cells to pool together to generate a simulated bulk sample. Default 500
#' @param min_chunkSize minimum number of cells required to construct the heter_cell_type
#' @param use_chunk a character indicating which cells to pool together for the heter_cell_type. Default='all' other options include 'random'
#'    When use_chunk = 'all', use all the cells belonging to the same patient for a given cell type to generate the certain cell type component in the simulated bulk;
#'    when use_chunk = 'random', randomly select 50-100% of the cells belonging to the same patient for a given cell type
#' @param export_cellUsage a logical variable determining whether to export cell names used to generate the simulated bulk. Default = F
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export

bulkSimulator_semi<-function(scExpr,scMeta,
                             colnames_of_cellType = NA,
                             colnames_of_sample = NA,
                             heter_cell_type = NA,
                             simulated_frac = NULL,
                             ncells_perSample=500,
                             min_chunkSize=10,
                             use_chunk = 'all',
                             export_cellUsage = F,
                             n.core = 1){
  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.na(colnames_of_sample)){
    stop('please provide column name that corresponds to sampleID in scMeta')
  }

  if(is.na(heter_cell_type)){
    stop('please provide the name of the cell type for which cross-sample heterogeneity should be maintained')
  }

  if(is.null(simulated_frac)){
    stop('please provide cell fraction matrix with pre-defined proportions for bulk simulation')
  }

  simulated_frac = as.matrix(simulated_frac)

  cell_type_labels=scMeta[,colnames_of_cellType]
  stopifnot(all(colnames(simulated_frac) %in% unique(cell_type_labels)))
  scMeta[,colnames_of_sample] = as.character(scMeta[,colnames_of_sample])

  if(!heter_cell_type %in% scMeta[, colnames_of_cellType]){
    stop('Please make sure the provided "heter_cell_type" is a cell type from scMeta.')
  }

  nCells_per_ct = table(scMeta[, colnames_of_cellType])
  if(unname(nCells_per_ct[heter_cell_type]) < min_chunkSize){
    warning("With the current min_chunkSize setting, the cell type-specific profiles for heter_cell_type will be identical across samples, as there are insufficient single cells available for subsetting ")
  }

  if(export_cellUsage){
    pboptions(type = "txt", style = 3, char = "=")

    sim_res = pblapply(seq(1,nrow(simulated_frac)),setCountSemi,simulated_frac,
                       scExpr,scMeta,
                       colnames_of_cellType,colnames_of_sample,
                       heter_cell_type,
                       ncells_perSample,
                       min_chunkSize,
                       use_chunk,
                       export_cellUsage,
                       cl = n.core)

    D = do.call(cbind, lapply(sim_res, function(x) x[[1]]))
    colnames(D)=paste0('simulated',1:ncol(D))

    CellNameList = lapply(sim_res, function(x) x[[2]])
    cellUsed = unlist(CellNameList)
    cell_usage = data.frame(cell = cellUsed,
                            bulk_id = paste0('simulated',rep(seq_along(CellNameList), lengths(CellNameList))))

    names(cell_type_labels) = rownames(scMeta)
    cell_usage$cell_type = cell_type_labels[match(cell_usage$cell,names(cell_type_labels))]

    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac,
                cell_usage = cell_usage))

  }else{
    pboptions(type = "txt", style = 3, char = "=")
    D=do.call(cbind,pblapply(seq(1,nrow(simulated_frac)),setCountSemi,simulated_frac,
                             scExpr,scMeta,
                             colnames_of_cellType,colnames_of_sample,
                             heter_cell_type,
                             ncells_perSample,
                             min_chunkSize,
                             use_chunk,
                             export_cellUsage,
                             cl = n.core))
    colnames(D)=paste0('simulated',1:ncol(D))

    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac))
  }
}

################### functions for heter simulation ###################
meta_shuffle<-function(scMeta){
  # Requires columns: "sampleID" and "cell_type"

  scMeta$sampleID = as.character(scMeta$sampleID)

  scMeta_shuffled=data.frame(matrix(NA,ncol = ncol(scMeta),nrow =0))

  cell_types=unique(scMeta$cell_type)

  for (i in 1:length(cell_types)){
    rows = which(scMeta$cell_type == cell_types[i])

    # replace sample labels
    meta_sub=scMeta[rows,]
    meta_sub=meta_sub[order(meta_sub$sampleID),]

    unique_samples=unique(meta_sub$sampleID)
    a = sample(unique_samples,length(unique_samples),replace = F)

    meta_sub$sampleID = rep(a,table(meta_sub$sampleID)%>% as.numeric())
    scMeta_shuffled=rbind(scMeta_shuffled,meta_sub)
  }
  return(scMeta_shuffled)
}

heter_aggregate<-function(ID,simulated_frac,scExpr,scMeta,colnames_of_cellType,
                          colnames_of_sample,
                          min_chunkSize = 5,
                          use_chunk = 'all',
                          export_cellUsage = F){

  simulated_frac_row = simulated_frac[ID,]

  # rename scMeta
  scMeta_renamed = data.frame(row.names = rownames(scMeta),
                              sampleID = scMeta[,colnames_of_sample],
                              cell_type = scMeta[,colnames_of_cellType])
  scMeta_shuffled = meta_shuffle(scMeta_renamed)

  sampleIDs = unique(scMeta_renamed$sampleID)
  selected_sample = sample(sampleIDs,1)

  cellType_component = names(simulated_frac_row)[simulated_frac_row > 0]
  scMeta_selected = subset(scMeta_shuffled, sampleID == selected_sample)
  scMeta_selected = scMeta_selected[scMeta_selected$cell_type %in% cellType_component,]

  scMeta_unselected = scMeta_shuffled[scMeta_shuffled$sampleID!=selected_sample,]

  p = rownames(scMeta_selected)

  chunk_size = scMeta_selected %>% group_by(cell_type) %>% summarise(n=n()) %>% deframe()

  missing_cellType = cellType_component[!cellType_component %in% names(chunk_size)]

  undersampled_cellType = names(chunk_size)[chunk_size < min_chunkSize]
  undersampled_cellType = intersect(undersampled_cellType,cellType_component)

  if(length(missing_cellType)>0){
    for(ct in missing_cellType){
      scMeta_candidate = scMeta_unselected[scMeta_unselected$cell_type == ct,]

      if(nrow(scMeta_candidate) < min_chunkSize){
        p = unique(c(p,rownames(scMeta_candidate)))
      }else{
        p = c(p,rownames(scMeta_candidate)[1:min_chunkSize])
      }
    }
  }

  if(length(undersampled_cellType) > 0 ){
    for(ct in undersampled_cellType){
      scMeta_candidate = scMeta_unselected[scMeta_unselected$cell_type == ct,]
      patch_n = min_chunkSize - chunk_size[ct] %>% as.numeric()

      if(nrow(scMeta_candidate) < patch_n){
        p = unique(c(p,rownames(scMeta_candidate)))
      }else{
        p = c(p,rownames(scMeta_candidate)[1:patch_n])
      }
    }
  }

  scMeta_final = scMeta_renamed[p,]
  group = list()
  for(ct in cellType_component){
    group[[ct]] <- rownames(scMeta_final)[scMeta_final$cell_type == ct]
  }

  if(use_chunk=='all'){
    X = scExpr[,p]
    C = do.call(cbind,lapply(group,function(x) Matrix::rowMeans(X[,x,drop=F])))
    E = C %*% matrix(simulated_frac_row[colnames(C)],ncol = 1)
  }else if(use_chunk=='random'){
    # randomly select 50%-100% cells in a chunk
    quick_select=function(cell_barcodes){
      percent=sample(seq(50,100),size = 1)
      selected=sample(cell_barcodes,size = ceiling(length(cell_barcodes)*percent*0.01))
      selected
    }
    group=lapply(group,quick_select)
    X = scExpr[,do.call(c,group)]
    C = do.call(cbind,lapply(group,function(x) Matrix::rowMeans(X[,x,drop=F])))
    E = C %*% matrix(simulated_frac_row[colnames(C)],ncol = 1)
  }
  E=apply(E,2,function(x)((x/sum(x))*1e+06))

  if(export_cellUsage){
    return(list(E = E,
                cellUsed = do.call(c,group)))
  }else{
    return(E)
  }
}

#' Heterogeneous bulk simulation
#' @description Generate bulk samples by aggregating single cells, with each cell-type component is constrained to cells originating from the same patient.
#'    This requires that sampleID is pre-defined in the scRNA dataset.
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns
#' @param min_chunkSize minimum number of cells required to construct a particular cell-type component in the simulated bulk, such as requiring at least 20 cells for B cells, at least 20 cells for T cells, and so forth.
#' @param use_chunk a character indicating which cells to pool together for the heter_cell_type. Default='all' other options include 'random'.
#'    When use_chunk = 'all', use all the cells belonging to the same patient for a given cell type to generate the certain cell type component in the simulated bulk;
#'    when use_chunk = 'random', randomly select 50-100% of the cells belonging to the same patient for a given cell type
#' @param export_cellUsage a logical variable determining whether to export cell names used to generate the simulated bulk. Default = F
#' @param n.core number of cores to use for parallel programming. Default = 1
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
bulkSimulator_heter <- function(scExpr, scMeta,
                                colnames_of_cellType = NA,
                                colnames_of_sample = NA,
                                simulated_frac = NULL,
                                min_chunkSize = 5,
                                use_chunk = 'all',
                                export_cellUsage = FALSE,
                                n.core = 1) {
  # Check input correspondence
  if (!all(colnames(scExpr) == rownames(scMeta))) {
    stop("colnames(scExpr) must match rownames(scMeta) (and in the same order).")
  }
  if (is.na(colnames_of_cellType) || !(colnames_of_cellType %in% colnames(scMeta))) {
    stop('Please provide a valid column name that corresponds to cell type in scMeta.')
  }
  if (is.na(colnames_of_sample) || !(colnames_of_sample %in% colnames(scMeta))) {
    stop('Please provide a valid column name that corresponds to sampleID in scMeta.')
  }
  if (is.null(simulated_frac)) {
    stop('Please provide cell fraction matrix with pre-defined proportions for bulk simulation.')
  }

  simulated_frac <- as.matrix(simulated_frac)
  scMeta[, colnames_of_sample] <- as.character(scMeta[, colnames_of_sample])
  scMeta[, colnames_of_cellType] <- as.character(scMeta[, colnames_of_cellType])

  nCells_per_ct = table(scMeta[, colnames_of_cellType])
  if(any(nCells_per_ct < min_chunkSize)){
    warning(
      paste0(
        "With the current min_chunkSize setting, the cell type-specific profiles for the following cell types will be identical across samples, as there are insufficient single cells available for subsetting: \n ",
        paste(names(nCells_per_ct)[nCells_per_ct < min_chunkSize], collapse = ", "),
        "\n Consider decreasing min_chunkSize or selecting a single-cell dataset with larger cell counts."
      )
    )
  }

  cell_type_labels <- scMeta[, colnames_of_cellType]
  if (!all(colnames(simulated_frac) %in% unique(cell_type_labels))) {
    stop('All columns of simulated_frac must be cell types found in scMeta.')
  }

  pboptions(type = "txt", style = 3, char = "=")
  if (export_cellUsage) {
    sim_res <- pblapply(seq_len(nrow(simulated_frac)), heter_aggregate,
                        simulated_frac,
                        scExpr,
                        scMeta,
                        colnames_of_cellType,
                        colnames_of_sample,
                        min_chunkSize,
                        use_chunk,
                        export_cellUsage,
                        cl = n.core)
    D <- do.call(cbind, lapply(sim_res, function(x) x[[1]]))
    colnames(D) <- paste0('simulated', seq_len(ncol(D)))
    CellNameList <- lapply(sim_res, function(x) x[[2]])
    cellUsed <- unlist(CellNameList)
    cell_usage <- data.frame(cell = cellUsed,
                             bulk_id = paste0('simulated', rep(seq_along(CellNameList), lengths(CellNameList))))

    names(cell_type_labels) <- rownames(scMeta)
    cell_usage$cell_type <- cell_type_labels[match(cell_usage$cell, names(cell_type_labels))]

    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac,
                cell_usage = cell_usage))
  } else {
    D <- do.call(cbind, pblapply(seq_len(nrow(simulated_frac)), heter_aggregate,
                                 simulated_frac,
                                 scExpr,
                                 scMeta,
                                 colnames_of_cellType,
                                 colnames_of_sample,
                                 min_chunkSize,
                                 use_chunk,
                                 cl = n.core))
    colnames(D) <- paste0('simulated', seq_len(ncol(D)))
    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac))
  }
}

############## functions for sampleID-free heter simulation ###########
#' sampleID-independent heterogeneous bulk simulation
#' @description  Generate bulk samples by aggregating cells from sub-clusters of each cellType. The simulated bulk sample comprises a single sub-cluster component
#'    for the heter_cell_type, with heter_cell_type typically recommended to be set as the malignant cell type. For the remaining cell type components, it includes random sub-cluster information
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_subcluster column name that corresponds to subcluster info in scMeta, where subcluster contains sub-clustering information for each cellType.
#'    Set to NA if subclustering information is not available; the function will then generate subclustering information automatically
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns
#' @param heter_cell_type name of the cell_type to maintain the highest level of heterogeneity. It is recommended to set this parameter to the name of the malignant cell types
#' @param export_cellUsage a logical variable determining whether to export cell names used to generate the simulated bulk. Default = F
#' @param dirichlet_cs_par a numeric value determine the dispersion level of the simulated fractions. With lower value indicating higher dispersion level. Default = 0.1
#' @param min.subcluster.size minimum size of a subcluster. This argument is required when colnames_of_subcluster is NA. Default = 20
#' @param max.num.cs maximum number of sub-clusters to aggregate for each cluster. If set to NA, this number will be equal to number of unique sub-clusters within each cell type
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
bulkSimulator_heter_sampleIDfree = function(scExpr,scMeta,
                                            colnames_of_cellType = NA,
                                            colnames_of_subcluster = NA,
                                            simulated_frac = NULL,
                                            heter_cell_type = NA,
                                            export_cellUsage = F,
                                            dirichlet_cs_par=0.1,
                                            min.subcluster.size = 20,
                                            max.num.cs = NA){

  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.na(heter_cell_type)){
    stop('please provide the name of the cell type for which cross-sample heterogeneity should be maintained')
  }

  if(is.null(simulated_frac)){
    stop('please provide cell fraction matrix with pre-defined proportions for bulk simulation')
  }

  cell_type_labels=scMeta[,colnames_of_cellType]
  stopifnot(all(colnames(simulated_frac) %in% unique(cell_type_labels)))

  if(!is.na(colnames_of_subcluster)){
    subcluster_labels = scMeta[,colnames_of_subcluster]
  }else{

    if (!requireNamespace("scran", quietly = TRUE)){
      stop('Please make sure you have scran installed to enable subclustering')
    }

    subcluster_labels = get_subcluster(scExpr,cell_type_labels, min.subcluster.size) %>% suppressWarnings()
  }

  C = build_ref_matrix(scExpr,subcluster_labels)
  meta = data.frame(cluster = cell_type_labels,
                    subcluster = subcluster_labels,
                    row.names = rownames(scMeta))

  simulate_dirichlet_cs_frac<-function(meta,ct_frac,dirichlet_cs_par){
    subcluster_table=meta %>% group_by(cluster,subcluster) %>% summarise(n=n()) %>% as.data.frame()
    rownames(subcluster_table)=subcluster_table$subcluster

    cell_types=unique(subcluster_table$cluster)
    cs_frac_all=matrix(NA,ncol =0,nrow = nrow(ct_frac))
    for(cell_type in cell_types){
      cs=subcluster_table$subcluster[subcluster_table$cluster==cell_type]
      cs_frac=matrix(NA,ncol = length(cs),nrow = nrow(ct_frac))
      colnames(cs_frac)=cs
      for (i in 1:nrow(ct_frac)){

        if(cell_type == heter_cell_type){
          num.cs = 1
        }else{
          if(is.na(max.num.cs)){
            num.cs=sample(1:length(cs),1,prob = rev(1:length(cs)))
          }else{
            max.num.cs.ct = min(max.num.cs,length(cs))
            num.cs=sample(1:max.num.cs.ct,1,prob = rev(1:max.num.cs.ct))
          }
        }
        cs.included=sample(cs,num.cs)
        cs_frac[i,cs.included]=gtools::rdirichlet(1,subcluster_table[cs.included,'n']*dirichlet_cs_par)[1,]
      }
      cs_frac[is.na(cs_frac)]=0
      cs_frac=sweep(cs_frac,1,ct_frac[,cell_type],'*')
      cs_frac_all=cbind(cs_frac_all,cs_frac)
    }
    return(cs_frac_all)
  }

  cs_frac = simulate_dirichlet_cs_frac(meta,simulated_frac,dirichlet_cs_par)

  D = C %*% t(cs_frac[,colnames(C)])
  D=apply(D,2,function(x)((x/sum(x))*1e+06))
  colnames(D)=paste0('simulated',1:ncol(D))

  if(export_cellUsage){
    rownames(cs_frac) = colnames(D)

    get_cell_used = function(bulk_id){
      non_zero_cs = colnames(cs_frac)[cs_frac[bulk_id,]!=0]
      out = data.frame(cell = rownames(meta)[meta$subcluster %in% non_zero_cs], bulk_id = bulk_id)
      return(out)
    }

    bulk_ids = colnames(D)

    cell_usage = do.call(rbind,lapply(bulk_ids,get_cell_used))
    rownames(cell_usage) = NULL

    names(cell_type_labels) = rownames(scMeta)
    cell_usage$cell_type = cell_type_labels[match(cell_usage$cell,names(cell_type_labels))]

    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac,
                cell_usage = cell_usage,
                meta = meta,
                cs_frac = cs_frac))

  }else{
    return(list(simulated_bulk = D,
                simulated_frac = simulated_frac))
  }

}


################ bulk simulation from other resources ##############
#' Bulk simulation using Generator() function from Favilaco et al.
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param nbulk number of simulated bulk samples
#' @param pool.size number of cells to aggregate for each simulated bulk sample
#' @param min.percentage minimum percentage of cellType fraction to generate in fraction simulation. Default = 1
#' @param max.percentage maximum percentage of cellType fraction to generate in fraction simulation. Default = 99
#' @param seed a single value. Default = 24
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
bulkSimulator_favilaco = function(scExpr, scMeta, colnames_of_cellType = NA, nbulk = 100, pool.size = 100,
                                  min.percentage = 1, max.percentage = 99, seed = 24){

  # this function incorporated the Generator() function from https://github.com/favilaco/deconv_benchmark/blob/master/helper_functions.R
  # the Generator() function expects the 'phenoData' argument to contain two columns named 'cellID' and 'cellType'
  # Furthermore, this function does not require a pre-defined proportion matrix, it will generate and provide the simulated fraction matrix autmatically

  warning('many functions used in the original favilaco methods are not longer supported in newer R version; this is a legacy version')

  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  phenoData <- scMeta %>% mutate_('cellType'=colnames_of_cellType)
  phenoData$cellID = rownames(phenoData)

  v = Generator(scExpr,phenoData,Num.mixtures = nbulk, pool.size, min.percentage, max.percentage, seed)
  simulated_frac = as.matrix(t(v$P))
  rownames(simulated_frac) = NULL

  D = as.matrix(v$T)
  D=apply(D,2,function(x)((x/sum(x))*1e+06))
  colnames(D)=paste0('simulated',1:ncol(D))

  return(list(simulated_bulk = D,
              simulated_frac = simulated_frac))
}


#' Bulk simulation using make_bulk_eset() function from immunedeconv package
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns
#' @param n_cells number of cells to aggregate for each simulated bulk sample
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
bulkSimulator_immunedeconv = function(scExpr, scMeta, colnames_of_cellType = NA , simulated_frac = NULL, n_cells = 500){

  require(immunedeconv, quietly = T)

  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.null(simulated_frac)){
    stop('please provide cell fraction matrix with pre-defined proportions for bulk simulation')
  }

  cell_type_labels=scMeta[,colnames_of_cellType]
  stopifnot(all(colnames(simulated_frac) %in% unique(cell_type_labels)))

  # immunedeconv::make_bulk_eset requires eset to have a `cell_type` column in `pData`.
  scMeta_renamed = data.frame(row.names = rownames(scMeta),
                              cell_type = scMeta[,colnames_of_cellType])

  fdata = data.frame(gene_symbol = rownames(scExpr))
  rownames(fdata) = rownames(scExpr)

  eset = Biobase::ExpressionSet(as.matrix(scExpr),phenoData = as(scMeta_renamed, "AnnotatedDataFrame"),featureData = as(fdata, "AnnotatedDataFrame"))
  new_eset <- immunedeconv::make_bulk_eset(eset, simulated_frac, n_cells = n_cells) %>% base::suppressMessages()
  D = exprs(new_eset)
  colnames(D)=paste0('simulated',1:ncol(D))

  return(list(simulated_bulk = D,
              simulated_frac = simulated_frac))
}


#' Bulk simulation using generateBulk_norep() function from SCDC package
#'
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta a dataframe that stores annotation info of each cells
#' @param colnames_of_cellType column name that corresponds to cellType annotation in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta
#' @param nbulk number of simulated bulk samples
#' @param disease indicate the health condition of subjects
#' @param ct.sub a subset of cell types that are selected to construct pseudo bulk samples. If NULL, then all cell types are used.
#' @param prop_mat manually input the cell-type proportion for pseudo bulk samples; with cellTypes in columns and samples in rows
#' @param samplewithRep logical, randomly sample single cells with replacement. Default is T
#'
#' @return a list containing the simulated bulk expression and its associated simulated fraction matrix
#' @export
bulkSimulator_SCDC = function(scExpr,scMeta,colnames_of_cellType = NA, colnames_of_sample = NA, nbulk = 10, disease = NULL, ct.sub = NULL, prop_mat = NULL,
                          samplewithRep = T){
  require(SCDC, quietly = T)
  scExpr = as.matrix(scExpr)
  eset = Biobase::ExpressionSet(assayData = scExpr,phenoData = new("AnnotatedDataFrame", data = scMeta))

  if(is.na(colnames_of_cellType)){
    stop('please provide column name that corresponds to cellType in scMeta')
  }

  if(is.na(colnames_of_sample)){
    stop('please provide column name that corresponds to sampleID in scMeta')
  }

  if(!is.null(prop_mat)){
    # generateBulk_norep() function requires the colnames of prop_mat to be ordered alphabetically
    prop_mat = prop_mat[,order(colnames(prop_mat))]

    # generateBulk_norep() also requires that all cellTypes in prop_mat must present in all patients
    q = table(scMeta[,colnames_of_sample],scMeta[,colnames_of_cellType])
    q = q[,colnames(prop_mat)]
    stopifnot(all(q>0))

    stopifnot(nrow(prop_mat)==nbulk)
  }

  if(is.null(ct.sub)){
    ct.sub = unique(scMeta[,colnames_of_cellType])
  }
  v = SCDC::generateBulk_norep(eset, ct.varname = colnames_of_cellType, sample = colnames_of_sample,
                               disease, ct.sub, prop_mat, nbulk, samplewithRep) %>% base::suppressMessages()

  D = v$pseudo_bulk
  D=apply(D,2,function(x)((x/sum(x))*1e+06))
  colnames(D)=paste0('simulated',1:ncol(D))

  frac = v$true_p
  rownames(frac) = NULL

  return(list(simulated_bulk = D,
              simulated_frac = frac))
}

##################### bulk simulation function ##################
#' Bulk simulation function
#' @description Generate a list of bulk simulation objects using user-selected simulation methods
#'
#' @param methods a character vector indicating which bulk simulation methods to use. Use list_bulkSimulator() to check for available method names.
#' @param scExpr single cell expression matrix used to simulate bulk data, with genes in rows and cells in columns
#' @param scMeta dataframe that stores annotation info of each cells, rownames of scMeta should be equal to colnames of scExpr
#' @param colnames_of_cellType column name that corresponds to cellType in scMeta
#' @param colnames_of_sample column name that corresponds to sampleID in scMeta. Required for 'semi', 'heter', and 'SCDC' methods.
#' @param simulated_frac a matrix with pre-defined fraction of different cell types, with samples in rows and cell_types in columns. This argument if required for 'homo', 'semi', 'heter', 'heter_sampleIDfree', 'immunedeconv' methods
#' @param heter_cell_type name of the cell_type to maintain the highest level of heterogeneity. It is recommended to set this parameter to the name of the malignant cell types.
#'     This argument is required for 'semi' and 'heter_sampleIDfree' methods
#' @param ncells_perSample number of cells to aggregate for each simulated bulk sample. This is an argument required for 'homo', 'semi', 'immunedeconv' and 'SCDC' methods
#' @param min_chunkSize minimum number of cells required to construct a particular cell-type component in the simulated bulk, such as requiring at least 20 cells for B cells, at least 20 cells for T cells, and so forth. This is an argument required for 'semi' and 'heter' methods
#' @param use_chunk a character indicating which cells to pool together for the a given cell_type. Default='all' other options include 'random'.
#'    When use_chunk = 'all', use all the cells belonging to the same patient for a given cell type to generate the certain cell type component in the simulated bulk;
#'    when use_chunk = 'random', randomly select 50-100% of the cells belonging to the same patient for a given cell type. This is an argument required for 'semi' and 'heter' methods
#' @param colnames_of_subcluster column name that corresponds to subcluster info in scMeta, where subcluster contains sub-clustering information for each cellType.  This is an argument required for 'heter_sampleIDfree' method only.
#'    Set to NA if subclustering information is not available; the function will then generate subclustering information automatically
#' @param export_cellUsage a logical variable determining whether to export cell names used to generate the simulated bulk. Default = F. This is an argument is only applicable to 'homo', 'semi', 'heter' and 'heter_sampleIDfree' methods
#' @param nbulk number of simulated bulk samples. This argument is only applicable to 'SCDC' methods, as they will need this argument to generate the fraction matrix. Default = 100
#' @param n.core number of cores to use for parallel programming. Default = 1. Parallel programming is applicable to 'homo', 'semi' and 'heter' methods.
#' @param ... additional arguments to be passed to the following functions: bulkSimulator_heter_sampleIDfree(),and bulkSimulator_SCDC()
#'
#' @return a list of simulated bulk objects
#' @export
#'
#' @examples
#' \dontrun{
#' # first create a fraction matrix for the simulated bulk
#' simulated_frac = fracSimulator_Beta(scMeta, n = 100,
#'                                     colnames_of_sample = 'sampleID',
#'                                     colnames_of_cellType = 'cell_type',
#'                                     fixed_cell_type = 'malignant')
#'
#' bulkSimulator(methods = c('homo', 'semi','heter'),
#'               scExpr = scExpr,
#'               scMeta = scMeta,
#'               colnames_of_cellType = 'cell_type',colnames_of_sample = 'sampleID',
#'
#'               # pre-defined simulated_frac will be used for 'homo', 'semi, 'heter' and 'immunedeconv' methods
#'               simulated_frac = simulated_frac,
#'
#'               # argument for semi/heter_sampleIDfree bulk simulation method
#'               heter_cell_type = 'malignant',
#'
#'               # general simulation parameters
#'               ncells_perSample = 500,
#'               min_chunkSize = 20,
#'               use_chunk = 'random',
#'
#'               n.core = 16)
#' }
bulkSimulator = function(methods,
                         scExpr,
                         scMeta,
                         colnames_of_cellType = NA,
                         colnames_of_sample = NA,
                         simulated_frac=NULL,
                         heter_cell_type = NA,
                         ncells_perSample=500,
                         min_chunkSize=5,
                         use_chunk = 'all',
                         colnames_of_subcluster = NA,
                         export_cellUsage = F,
                         nbulk = 100,
                         n.core = 1,
                         ...){

  arg_list = list(...)
  heter_sampleIDfree_args = arg_list[names(arg_list) %in% c('dirichlet_cs_par','min.subcluster.size','max.num.cs')]
  SCDC_args = arg_list[names(arg_list) %in% c('disease','ct.sub','prop_mat','samplewithRep')]

  if('favilaco' %in% methods){
    stop('favilaco bulk simulation methods is no longer supported since deconvBenchmarking v0.1.1')
  }

  supported_methods <- c('homo', 'semi', 'heter', 'heter_sampleIDfree', 'immunedeconv', 'SCDC')
  if (any(!methods %in% supported_methods)) {
    stop(
      "Some bulkSimulation methods are not supported. Please provide only valid names. Supported methods are: ",
      paste(supported_methods, collapse = ", ")
    )
  }

  if('heter_sampleIDfree' %in% methods){

    if (!requireNamespace("scran", quietly = TRUE)){
      stop('Please make sure you have scran installed to run heter_sampleIDfree bulk simulation')
    }
  }

  if('immunedeconv' %in% methods){
    if (!requireNamespace("immunedeconv", quietly = TRUE)){
      stop('Please make sure you have immunedeconv installed to run immunedeconv bulk simulation')
    }
  }

  if('SCDC' %in% methods){
    if (!requireNamespace("SCDC", quietly = TRUE)){
      stop('Please make sure you have SCDC installed to run SCDC bulk simulation')
    }
  }

  if (any(c('heter', 'semi') %in% methods)) {
    if (is.na(heter_cell_type)) {
      stop('Please provide "heter_cell_type" when using semi/heter simulation.')
    } else {
      if (!heter_cell_type %in% unique(scMeta[,colnames_of_cellType])) {
        stop('Please make sure the provided "heter_cell_type" is a cell type from scMeta.')
      }
    }
  }

  bulk_list = list()
  l = list_bulkSimulator(show_description=T)

  for(method in methods){
    description = l$description[l$method == method]
    message(description)

    switch(method,
           homo = {
             result = bulkSimulator_homo(scExpr,scMeta,
                                         colnames_of_cellType,
                                         simulated_frac,
                                         ncells_perSample,
                                         export_cellUsage,
                                         n.core)
           },
           semi = {
             result = bulkSimulator_semi(scExpr,scMeta,
                                         colnames_of_cellType,
                                         colnames_of_sample,
                                         heter_cell_type,
                                         simulated_frac,
                                         ncells_perSample,
                                         min_chunkSize,
                                         use_chunk,
                                         export_cellUsage,
                                         n.core)
           },
           heter = {
             result = bulkSimulator_heter(scExpr,scMeta,
                                          colnames_of_cellType,
                                          colnames_of_sample,
                                          simulated_frac,
                                          min_chunkSize,
                                          use_chunk,
                                          export_cellUsage,
                                          n.core)
           },
           heter_sampleIDfree = {
             result = do.call(bulkSimulator_heter_sampleIDfree,c(list(scExpr = scExpr,
                                                                      scMeta = scMeta,
                                                                      colnames_of_cellType = colnames_of_cellType,
                                                                      colnames_of_subcluster = colnames_of_subcluster,
                                                                      simulated_frac = simulated_frac,
                                                                      heter_cell_type = heter_cell_type,
                                                                      export_cellUsage = export_cellUsage),
                                                                 heter_sampleIDfree_args))

           },
           immunedeconv = {
             result = bulkSimulator_immunedeconv(scExpr, scMeta,
                                                 colnames_of_cellType,
                                                 simulated_frac,
                                                 n_cells = ncells_perSample)
           },
           SCDC = {

             result = do.call(bulkSimulator_SCDC,c(list(scExpr = scExpr,
                                                        scMeta = scMeta,
                                                        colnames_of_cellType = colnames_of_cellType,
                                                        colnames_of_sample = colnames_of_sample,
                                                        nbulk = nbulk),
                                                   SCDC_args))
           },
           {
             warning(paste0("Invalid method specified: ",method, ", please use list_bulkSimulator() to check for available methods"))
             result <- NULL
           })
    bulk_list <- c(bulk_list, list(result))
  }
  names(bulk_list) = methods

  return(bulk_list)
}


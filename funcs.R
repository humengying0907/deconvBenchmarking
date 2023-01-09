#################### helper functions #####################################
commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}

scale01=function(x){x=x-min(x);x/max(x)}

build_ref_matrix<-function(Expr,cell_type_labels){
  stopifnot(ncol(Expr)==length(cell_type_labels))
  group = list()
  for(i in unique(cell_type_labels)){ 
    group[[i]] <- which(cell_type_labels %in% i)
  }
  C = do.call(cbind, lapply(group,function(x) Matrix::rowMeans(Expr[,x,drop=F])))
  C
}

downsampling<-function(scMeta,colnames_of_sample,colnames_of_cellType,scale=0.2){
  scMeta_summary=scMeta %>% group_by_(colnames_of_sample, colnames_of_cellType) %>%
    summarise(n=n()) %>%
    mutate(freq=n/sum(n)) %>%group_by_(colnames_of_cellType) %>% summarise(sum=sum(n))
  expected_n=floor(nrow(scMeta)*scale)
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

quick_cv<-function(v){
  sd(v)/mean(v)
}

make_unique<-function(v){
  return(v[!duplicated(v)])
}

get_cell_map<-function(performance_obj){
  build_df=function(deconv_name,deconv_results){
    data.frame(method=deconv_name,method_cell_type=rownames(deconv_results[[deconv_name]]))
  }
  deconv_names=names(performance_obj[[1]]$homo_deconvResults)
  cell_map=do.call(rbind,lapply(deconv_names,build_df,performance_obj[[1]]$homo_deconvResults)) 
  cell_map$ID=seq(1,nrow(cell_map))
  cell_map
}

get_mapping_list<-function(cell_map,simulated_frac){
  mapping=list()
  for(i in colnames(simulated_frac)){
    mapping[[i]]=which(cell_map$method_cell_type==i)
  }
  names(mapping)=colnames(simulated_frac)
  
  linseed_rowID=cell_map$ID[cell_map$method=='linseed']
  for(i in 1:length(mapping)){
    mapping[[i]]=c(mapping[[i]],linseed_rowID)
  }
  CAMfree_rowID=cell_map$ID[cell_map$method=='CAMfree']
  for(i in 1:length(mapping)){
    mapping[[i]]=c(mapping[[i]],CAMfree_rowID)
  }
  mapping
}


#################### functions used for cell fraction simulation ####################
# simulated from beta distribution of cell fractions in the single cell profile
simulate_frac<-function(scMeta,n,colnames_of_sample,colnames_of_cellType,fixed_cell_type=NULL){
  # n: number of simulated sampes
  # fixed_cell_type: rbeta for fixed_cell_type will not be scaled, fixed_cell_type is usually the most abundant cell types
  ct_table=scMeta %>% group_by_(colnames_of_sample, colnames_of_cellType) %>%
    summarise(n=n()) %>%
    mutate(freq=n/sum(n))
  
  if(is.null(fixed_cell_type)){
    fixed_cell_type=names(table(scMeta$cell_type))[which.max(table(scMeta$cell_type))]
  }
  
  # for fixed cell type, simulate cell fraction based on beta distribution
  x=ct_table[,'freq'][ct_table[,colnames_of_cellType]==fixed_cell_type]
  fit_x <- fitdist(x, "beta")
  fixed_frac=rbeta(n,fit_x$estimate[1],fit_x$estimate[2])
  
  quick_hist=function(){
    par(mfrow=c(1,2))
    hist(x,main = 'original frac for fixed cell type')
    hist(fixed_frac,main = 'simulated frac for fixed cell type')
  }
  quick_hist() # make sure that for fixed cell type, the distribution is similar
  
  # for other cell types, simulate cell fraction based on beta distribution, and scaled to the remaining frac
  quick_rbeta<-function(cell_type){
    x=ct_table[,'freq'][ct_table[,colnames_of_cellType]==cell_type]  
    if(1 %in% x){
      x[x==1]=x[x==1]-0.01
    }
    fit_x <- fitdist(x, "beta") # fitdist can fail when: https://stackoverflow.com/questions/44507568/error-when-fitting-a-beta-distribution-the-function-mle-failed-to-estimate-the  
    rbeta(n,fit_x$estimate[1],fit_x$estimate[2])
  }
  
  all_cell_types=unique(ct_table[,colnames_of_cellType])
  other_cell_types=all_cell_types[all_cell_types!=fixed_cell_type]
  flexible_frac=do.call(cbind,lapply(other_cell_types,quick_rbeta))
  
  remaining_frac=1-fixed_frac
  scaling_factor=remaining_frac/rowSums(flexible_frac)
  flexible_frac=sweep(flexible_frac,1,scaling_factor,'*')
  
  simlated_frac=cbind(fixed_frac,flexible_frac)
  colnames(simlated_frac)=c(fixed_cell_type,other_cell_types)
  simlated_frac
}

################### functions for heter simulation (using chunks) ###################
meta_shuffle<-function(scMeta,
                       colnames_of_cellType,
                       colnames_of_sample){
  scMeta_shuffled=data.frame(matrix(NA,ncol = ncol(scMeta),nrow =0))
  cell_types=unique(scMeta[,colnames_of_cellType])
  for (i in 1:length(cell_types)){
    rows=which(scMeta[,colnames_of_cellType]==cell_types[i])
    # replace sample labels
    meta_sub=scMeta[rows,]
    meta_sub=meta_sub[order(meta_sub[,colnames_of_sample]),]
    
    unique_samples=unique(meta_sub[,colnames_of_sample])
    a=sample(unique_samples,length(unique_samples),replace = F)
    
    meta_sub[,colnames_of_sample]=rep(a,table(meta_sub[,colnames_of_sample])%>% as.numeric())
    scMeta_shuffled=rbind(scMeta_shuffled,meta_sub)
  }
  return(scMeta_shuffled)
}

heter_aggregate<-function(ID,frac_table,scExpr,scMeta,colnames_of_cellType,
                          colnames_of_sample,
                          chunk_size_threshold=5,
                          use_chunk='all',
                          scale_to_million=T){
  
  # frac_table: simulated fraction with samples in rows and cell types in columns
  # scExpr: genes in rows and cell_barcode in columns
  # scMeta: dataframe that stores scMeta info of each cells, rownames of scMeta should be equal to colnames of scExpr
  # chunk_size_threshold: the minimum number of cells to aggregate for every cell type
  # use_chunk: whether or not use all chunk to create reference profile, default='all'; other options include 'random'
  
  stopifnot(all.equal(rownames(scMeta),colnames(scExpr)))
  
  simulated_frac=frac_table[ID,]
  scMeta=meta_shuffle(scMeta,colnames_of_cellType,colnames_of_sample)
  
  sampleIDs=unique(scMeta[,colnames_of_sample])
  patient_barcode=sample(sampleIDs,1)
  
  p=rownames(scMeta)[scMeta[,colnames_of_sample]==patient_barcode]
  
  # borrow a chunk if scMeta[p,] does not have enough cell types
  unique_cellType = unique(scMeta[,colnames_of_cellType])
  scMeta_sub=scMeta[p,]
  missing_cellType=unique_cellType[!unique_cellType %in% unique(scMeta_sub[,colnames_of_cellType])]
  
  if(length(missing_cellType)>0){
    for (i in 1:length(missing_cellType)){
      scMeta_borrow = scMeta[scMeta[,colnames_of_cellType]==missing_cellType[i],]
      unique_patient = unique(scMeta_borrow[,colnames_of_sample]) 
      borrow_patient = sample(unique_patient,1)
      # change borrow_patient id to patient barcode
      scMeta[scMeta[,colnames_of_cellType]==missing_cellType[i] & scMeta[,colnames_of_sample]==borrow_patient,colnames_of_sample]=patient_barcode 
    }
  }
  
  p=rownames(scMeta)[scMeta[,colnames_of_sample]==patient_barcode]
  
  # combine with random chunks if scMeta[p,cell_type] has less than chunk_size_threshold cells
  scMeta_sub=scMeta[p,]
  chunk_size=scMeta_sub %>% group_by_(colnames_of_cellType) %>% summarise(n=n()) %>% deframe() 
  
  undersampled_cellType=names(chunk_size)[chunk_size<chunk_size_threshold]
  if(length(undersampled_cellType)>0){
    for (i in 1:length(undersampled_cellType)){
      scMeta_borrow = scMeta[scMeta[,colnames_of_cellType]==undersampled_cellType[i],]
      unique_patient = unique(scMeta_borrow[,colnames_of_sample]) 
      
      n_after_borrow = chunk_size[undersampled_cellType[i]] %>% as.numeric()
      borrow_query=0
      while (n_after_borrow < chunk_size_threshold) {
        borrow_patient = sample(unique_patient,1)
        scMeta_borrow_patient_sub=scMeta_borrow[scMeta_borrow[,colnames_of_sample]==borrow_patient,]
        n_after_borrow=chunk_size[undersampled_cellType[i]] %>% as.numeric()+nrow(scMeta_borrow_patient_sub) # number of chunk size after borrow
        borrow_query=borrow_query+1
        if(borrow_query>10){
          stop('chunk_size_threshold too large: no enough cells to aggreate; use smaller chunk size threshold')
        }
      }
      scMeta[scMeta[,colnames_of_cellType]==undersampled_cellType[i] & scMeta[,colnames_of_sample]==borrow_patient,colnames_of_sample]=patient_barcode 
    }
  }
  
  # final cells to use
  p=rownames(scMeta)[scMeta[,colnames_of_sample]==patient_barcode] 
  scMeta_sub=scMeta[p,]
  chunk_size=scMeta_sub %>% group_by_(colnames_of_cellType) %>% summarise(n=n()) %>% deframe()
  # print(chunk_size)
  
  # generate simulated expression by matrix multiplication
  cellType <- scMeta[p,colnames_of_cellType]
  stopifnot(length(unique(cellType))==length(unique_cellType))
  group = list()
  for(i in unique(cellType)){ 
    group[[i]] <- p[which(cellType %in% i)]
  }
  if(use_chunk=='all'){
    # use all chunk information to create reference profile
    X = scExpr[,p]
    C = lapply(group,function(x) Matrix::rowMeans(X[,x,drop=F])) 
    C = do.call(cbind, C)
    E = C %*% matrix(simulated_frac[colnames(C)],ncol = 1)
  }else if(use_chunk=='random'){
    # randomly select 50%-100% cells of each chunk to create reference profile
    quick_select=function(cell_barcodes){
      percent=sample(seq(50,100),size = 1)
      selected=sample(cell_barcodes,size = ceiling(length(cell_barcodes)*percent*0.01))
      selected
    }
    group=lapply(group,quick_select)
    X = scExpr[,do.call(c,group)]
    C = lapply(group,function(x) Matrix::rowMeans(X[,x,drop=F])) 
    C = do.call(cbind, C)
    E = C %*% matrix(simulated_frac[colnames(C)],ncol = 1)
  }
  if (scale_to_million==T){
    E=apply(E,2,function(x)((x/sum(x))*1e+06)) 
  }
  return(E)
}


create_heterSimulation<-function(frac_table,scExpr,scMeta,colnames_of_cellType,
                                 colnames_of_sample,
                                 chunk_size_threshold=5,
                                 use_chunk='all',
                                 scale_to_million=T){
  D=do.call(cbind,lapply(seq(1,nrow(frac_table)),heter_aggregate,
                         frac_table,
                         scExpr,scMeta,
                         colnames_of_cellType,colnames_of_sample,
                         chunk_size_threshold,use_chunk,scale_to_million))
  colnames(D)=paste0('simulated',1:ncol(D))
  return(D)
}

######################## functions for homo simulation ############################
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
  # print(length(sc))
  return(sc)
}

setBulks<-function(ID,simulated_frac,CellNameList,scExpr,cell_type_labels){
  named_labels=cell_type_labels
  names(named_labels)=colnames(scExpr)
  
  simulated_frac_row=simulated_frac[ID,]
  p=CellNameList[[ID]]
  
  X=scExpr[,p]
  cell_type_labels_sub=named_labels[p] %>% unname()
  
  C=build_ref_matrix(X,cell_type_labels_sub)
  E=C %*% matrix(simulated_frac_row[colnames(C)],ncol = 1)
  E=apply(E,2,function(x)((x/sum(x))*1e+06)) 
  return(E)
}


create_homoSimulation<-function(scExpr,scMeta,colnames_of_cellType,simulated_frac,ncells_perSample=500){
  stopifnot(all.equal(colnames(scExpr),rownames(scMeta)))
  cell_type_labels=scMeta[,colnames_of_cellType]
  cellNames=colnames(scExpr)
  cellSet=list()
  for (ts in levels(factor(cell_type_labels))) {
    cellSet[[ts]] <- cellNames[cell_type_labels == ts]
  }
  CellNameList <- apply(simulated_frac,1,setCells,cellSet,ncells_perSample) %>% t() 
  D=do.call(cbind,lapply(seq(1,nrow(simulated_frac)),
                                 setBulks,simulated_frac,CellNameList,scExpr,cell_type_labels))
  colnames(D)=paste0('simulated',1:ncol(D))
  return(D)
}


######################### functions for semiheter simulation #####################
setCountSemi <- function(ID,frac_table,scExpr,scMeta,colnames_of_cellType,colnames_of_sample,fixed_cell_type,
                         ncells_perSample=500,chunk_size_threshold_for_fixed_cell=10){
  
  simulated_frac=frac_table[ID,]
  
  # fixed cell type is restricted to one patient 
  fixed_metaSub=scMeta[scMeta[,colnames_of_cellType]==fixed_cell_type,]
  sampleIDs=unique(fixed_metaSub[,colnames_of_sample])
  
  fixed_cells='initial_cell'
  query=0
  while(length(fixed_cells)<chunk_size_threshold_for_fixed_cell){
    fixed_sample=sample(sampleIDs,1)
    fixed_cells=rownames(fixed_metaSub)[fixed_metaSub[,colnames_of_sample]==fixed_sample]
    query=query+1
    if(query>10){
      stop('chunk_size_threshold_for_fixed_cell too large: no enough cells to aggreate; use smaller chunk size threshold')
    }
  }
  
  # non-malignant cells drawn from all patients
  cellNames=rownames(scMeta)
  cell_type_labels=scMeta[,colnames_of_cellType]
  cellSet=list()
  for (ts in levels(factor(cell_type_labels))) {
    cellSet[[ts]] <- cellNames[cell_type_labels == ts]
  }
  
  cellSet_remaining=cellSet[names(cellSet)!=fixed_cell_type]
  simulated_frac_remaining=simulated_frac[names(simulated_frac)!=fixed_cell_type]
  
  x=simulated_frac_remaining*ncells_perSample
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
  
  scMeta_sub=scMeta[sc,]
  chunk_size=scMeta_sub %>% group_by_(colnames_of_cellType) %>% summarise(n=n()) %>% deframe()
  
  # generate simulated expression by matrix multiplication
  cellType <- scMeta[sc,colnames_of_cellType]
  group = list()
  for(i in unique(cellType)){ 
    group[[i]] <- which(cellType %in% i)
  }
  
  X = scExpr[,sc]
  C = lapply(group,function(x) Matrix::rowMeans(X[,x,drop=F])) 
  C = do.call(cbind, C)
  
  E=C %*% matrix(simulated_frac[colnames(C)],ncol = 1)
  E=apply(E,2,function(x)((x/sum(x))*1e+06)) 
  
  return(E)
  
}

create_semiheterSimulation<-function(frac_table,scExpr,scMeta,colnames_of_cellType,
                                     colnames_of_sample,fixed_cell_type,
                                     ncells_perSample=500,chunk_size_threshold_for_fixed_cell=10){
  D=do.call(cbind,lapply(seq(1,nrow(frac_table)),setCountSemi,frac_table,scExpr,scMeta,
                         colnames_of_cellType,colnames_of_sample,fixed_cell_type,
                         ncells_perSample,chunk_size_threshold_for_fixed_cell))
  colnames(D)=paste0('simulated',1:ncol(D))
  return(D)
}

################## functions for DE analysis ##########################
# marker.fc() function comes from https://github.com/favilaco/deconv_benchmark
marker.fc <- function(fit2, cont.matrix,log2.threshold = 1, output_name = "markers"){
  
  topTable_RESULTS = limma::topTable(fit2, coef = 1:ncol(cont.matrix), number = Inf, adjust.method = "BH", p.value = 0.05, lfc = log2.threshold)
  AveExpr_pval <- topTable_RESULTS[,(ncol(topTable_RESULTS)-3):ncol(topTable_RESULTS)]
  topTable_RESULTS <- topTable_RESULTS[,1:(ncol(topTable_RESULTS)-4)]
  
  if(length(grep("ERCC-",topTable_RESULTS$gene)) > 0){ topTable_RESULTS <- topTable_RESULTS[-grep("ERCC-",topTable_RESULTS$gene),] }
  
  markers <- apply(topTable_RESULTS,1,function(x){
    temp = sort(x)
    ((temp[ncol(topTable_RESULTS)] - temp[ncol(topTable_RESULTS)-1]) >= log2.threshold) | (abs(temp[1] - temp[2]) >= log2.threshold)
    
  })
  
  topTable_RESULTS = topTable_RESULTS[markers,]
  
  markers <- cbind.data.frame(rownames(topTable_RESULTS),
                              t(apply(topTable_RESULTS, 1, function(x){
                                temp = max(x)
                                if(temp < log2.threshold){
                                  temp = c(min(x),colnames(topTable_RESULTS)[which.min(x)])
                                } else {
                                  temp = c(max(x),colnames(topTable_RESULTS)[which.max(x)])
                                } 
                                temp
                              })))
  
  colnames(markers) <- c("gene","log2FC","CT")
  markers$log2FC = as.numeric(as.character(markers$log2FC))
  markers <- markers %>% dplyr::arrange(CT,desc(log2FC)) 
  
  markers$AveExpr <- AveExpr_pval$AveExpr[match(markers$gene,rownames(AveExpr_pval))]
  markers$gene <- as.character(markers$gene)
  markers$CT <- as.character(markers$CT)
  
  #write.table(markers, file = output_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  return(markers)
  
}

one_against_all_statistics<-function(scExpr,cell_type_labels,hv_genes){
  # scExpr has genes in rows and samples in columns, and should be un-log transformed
  # cell_type_labels correspond to columns in scExpr
  stopifnot(ncol(scExpr)==length(cell_type_labels))
  stopifnot(max(scExpr)>100)
  annotation=factor(cell_type_labels)
  design <- model.matrix(~0+annotation)
  colnames(design) <- unlist(lapply(strsplit(colnames(design),"annotation"), function(x) x[2]))
  cont.matrix <- matrix((-1/ncol(design)),nrow=ncol(design),ncol=ncol(design))
  colnames(cont.matrix) <- colnames(design)
  diag(cont.matrix) <- (ncol(design)-1)/ncol(design)
  
  hv_genes=hv_genes[hv_genes%in% rownames(scExpr)]
  fit <- limma::lmFit(log2(scExpr[hv_genes,]+1), design)
  fit2 <- limma::contrasts.fit(fit, cont.matrix)
  fit2 <- limma::eBayes(fit2, trend=TRUE)
  
  return(list(fit2,cont.matrix)) 
}

pairwiseTTest_statistics<-function(scExpr,cell_type_labels,cell_state_labels,hv_genes){
  # scExpr has genes in rows and samples in columns, and should be un-log transformed
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

build_DE_standard<-function(DE_statistics,log2FC=2,minimum_n=15,maximum_n=50){
  fit2=DE_statistics[[1]]
  cont.matrix=DE_statistics[[2]]
  markers_loose=marker.fc(fit2,cont.matrix,log2.threshold = 1) # set up loose threshold in case that there's not enough genes for a cell type
  markers_loose=markers_loose[markers_loose$log2FC>0,]
  
  markers=marker.fc(fit2,cont.matrix,log2.threshold = log2FC)
  markers=markers[markers$log2FC>0,]
  
  DE_list<-list() 
  sig_genes=NULL
  DE_cellTypes=unique(markers$CT)
  
  for (i in 1:length(DE_cellTypes)){
    first_try=rownames(markers)[markers$CT==DE_cellTypes[i]]
    second_try=rownames(markers_loose)[markers_loose$CT==DE_cellTypes[i]]
    if(length(first_try)>=minimum_n){
      if(length(first_try)>maximum_n){
        message(paste('For',DE_cellTypes[i],':',length(first_try),'genes passed log2FC threshold, will pick top',maximum_n))
        DE_list[[i]]=first_try[1:maximum_n]
      }else{
        DE_list[[i]]=first_try
      }
    }else if (length(second_try)>=minimum_n){
      DE_list[[i]]=second_try[1:minimum_n]
      message(paste('For',DE_cellTypes[i],': not enough number of genes passing log2FC threshod, choose first',minimum_n,'genes that passed log2FC=1 threshold'))
    }else{
      DE_list[[i]]=second_try
      warning(paste('For',DE_cellTypes[i],': not enough number of genes passing log2FC=1 threshold, choose all genes passing log2FC=1 threshold'))
    }
    
    sig_genes=append(sig_genes,DE_list[[i]])
    names(DE_list)[i]=DE_cellTypes[i]
  }
  return(list(DE_list=DE_list,sig_genes=sig_genes))
}

build_DE_pairwiseTTest<-function(DE_statistics,log2FC=2,minimum_n=15,maximum_n=50){
  DE_list<-list() 
  sig_genes=NULL
  DE_cellTypes=names(DE_statistics)
  
  for (i in 1:length(DE_cellTypes)){
    sub=DE_statistics[[i]] %>% dplyr::arrange(desc(-pval.up.min)) 
    first_try=rownames(sub)[sub$min.lfc>2]
    second_try=rownames(sub)[sub$min.lfc>1]
    
    if(length(first_try)>=minimum_n){
      if(length(first_try)>maximum_n){
        message(paste('For',DE_cellTypes[i],':',length(first_try),'genes passed log2FC threshold, will pick top',maximum_n,'with lowest p-value'))
        DE_list[[i]]=first_try[1:maximum_n]
      }else{
        DE_list[[i]]=first_try
      }
    }
    else if (length(second_try)>=minimum_n){
      DE_list[[i]]=second_try[1:minimum_n]
      message(paste('For',DE_cellTypes[i],': not enough number of genes passing log2FC threshod, choose first',minimum_n,'genes that passed log2FC=1 threshold'))
    }else{
      DE_list[[i]]=second_try
      warning(paste('For',DE_cellTypes[i],': not enough number of genes passing log2FC=1 threshold, choose all genes passing log2FC=1 threshold'))
    }
    
    sig_genes=append(sig_genes,DE_list[[i]])
    names(DE_list)[i]=DE_cellTypes[i]
    
  }
  return(list(DE_list=DE_list,sig_genes=sig_genes))
}

################## function used for training/testing ############################
create_train_test_splitting<-function(scMeta,colnames_of_sample,colnames_of_cellType){
  # make sure that in training cells contain around half of each cell types
  # failing_criteria controls that number of cells is within 40%-60% percent of the original meta
  scMeta_summary=scMeta %>%   group_by_(colnames_of_sample, colnames_of_cellType) %>%
    summarise(n=n()) %>%
    mutate(freq=n/sum(n)) %>%group_by_(colnames_of_cellType) %>% summarise(sum=sum(n))
  
  repeat{
    training_cells=sample(rownames(scMeta),size = floor(nrow(scMeta)/2),replace = F)
    training_meta=scMeta[training_cells,]
    
    trainMeta_summary=training_meta %>% group_by_(colnames_of_sample, colnames_of_cellType) %>%
      summarise(n=n()) %>%
      mutate(freq=n/sum(n)) %>%group_by_(colnames_of_cellType) %>% summarise(sum=sum(n))
    
    failing_criteria=sum((trainMeta_summary$sum <=  scMeta_summary$sum*0.4) + (trainMeta_summary$sum >=  scMeta_summary$sum*0.6 ))
    
    if(failing_criteria==0){
      break
    }
  }
  
  testing_cells=rownames(scMeta)[!rownames(scMeta)%in%training_cells]
  
  splitting_list=list()
  splitting_list[['training_cells']]=training_cells
  splitting_list[['testing_cells']]=testing_cells
  splitting_list
  
}

create_simulation_fold <- function(scExpr,scMeta,colnames_of_sample,colnames_of_cellType,
                                     simulated_frac, # simulated cell fractions with samples in rows and cell types in columns
                                     max.spec_cutoff_for_DE=0.3, # parameters for hv genes 
                                     chunk_size_threshold=5,scale_to_million=T, # heter-simulation parameters
                                     ncells_perSample=500, # homo-simulation parameter
                                     fixed_cell_type,chunk_size_threshold_for_fixed_cell=10, # semiheter-simulation parameters
                                     log2FC=2,minimum_n=15,maximum_n=50, colnames_of_cellState=NULL,# marker gene list parameters (maximum_n is set higher in case of downsampling)
                                     create_autogeneS_input=F,autogeneS_input_file_name=NULL, max.spec_cutoff_for_autogeneS=0.5, # parametner for hv genes# create input for autogeneS
                                     create_cibersortx_input=F,cibersort_input_file_name=NULL,cibersort_downsample=F,cibersort_downsample_scale=0.2){
  stopifnot(nrow(scMeta)==ncol(scExpr))
  if(is.null(colnames_of_cellState)){
    colnames_of_cellState=colnames_of_cellType
  }
  
  simulation_fold=new.env()
  
  message('>>> split single cell object into training and testing cells')
  splitting=list()
  cell_splitting=create_train_test_splitting(scMeta ,colnames_of_sample,colnames_of_cellType)
  training_cells=cell_splitting$training_cells
  testing_cells=cell_splitting$testing_cells
  
  splitting$training_cells=training_cells
  splitting$testing_cells=testing_cells
  simulation_fold$splitting=splitting
  
  scMeta_train=scMeta[training_cells,]
  scMeta_test=scMeta[testing_cells,]
  
  scExpr_train=scExpr[,training_cells]
  scExpr_test=scExpr[,testing_cells]
  
  cell_type_labels_train=scMeta_train[,colnames_of_cellType]
  cell_type_labels_test=scMeta_test[,colnames_of_cellType]
  cell_state_labels_train=scMeta_train[,colnames_of_cellState]
  
  stopifnot(all.equal(rownames(scMeta_train),colnames(scExpr_train)))
  stopifnot(all.equal(rownames(scMeta_test),colnames(scExpr_test)))
  
  message('>>> find hv genes from scExpr_train')
  sc.stat <- BayesPrism::plot.scRNA.outlier(
    input=scExpr_train %>% t(), 
    cell.type.labels=cell_type_labels_train,
    species="hs", 
    return.raw=TRUE #return the data used for plotting
  )
  # max.spec is the maximum specificity of each gene
  hv_genes_for_DE=rownames(sc.stat)[sc.stat$max.spec>=max.spec_cutoff_for_DE]
  print(paste('use',length(hv_genes_for_DE),'genes for DE analysis'))

  message('>>> running DE analysis on scExpr_train')
  limma_statistics=one_against_all_statistics(scExpr_train,cell_type_labels_train,hv_genes_for_DE)
  scran_statistics=pairwiseTTest_statistics(scExpr_train,cell_type_labels_train,cell_state_labels_train,hv_genes_for_DE)
  DE_statistics=new.env()
  DE_statistics$DE_statistics_limma=limma_statistics
  DE_statistics$DE_statistics_scran=scran_statistics
  simulation_fold$DE_statistics=DE_statistics
  
  message('>>> summarizing DE statistics')
  DE_list_bundles=list()
  print('build DE list based on limma statistics')
  limma_results=build_DE_standard(limma_statistics,log2FC,minimum_n,maximum_n)
  print('build DE list based on scran statistics')
  scran_results=build_DE_pairwiseTTest(scran_statistics,log2FC,minimum_n,maximum_n)
  DE_list_bundles[['limmaMarkers']]=limma_results$DE_list
  DE_list_bundles[['scranMarkers']]=scran_results$DE_list
  simulation_fold$DE_list_bundles=DE_list_bundles
  
  message('>>> build reference Matrix')
  C=build_ref_matrix(scExpr_train,cell_type_labels_train)
  if(create_autogeneS_input==T){
    hv_genes_for_autogeneS=rownames(sc.stat)[sc.stat$max.spec>=max.spec_cutoff_for_autogeneS]
    print(paste('use',length(hv_genes_for_autogeneS),'genes for autogeneS input'))
    stopifnot(length(hv_genes_for_autogeneS)>2000) 
    write.table(C[hv_genes_for_autogeneS,],file = autogeneS_input_file_name,sep = ',') # input for autogeneS
  }
  ref_list=list()
  ref_list[['limmaRef']]=C[limma_results$sig_genes,]
  ref_list[['scranRef']]=C[scran_results$sig_genes,]
  simulation_fold$ref_list=ref_list
  
  if(create_cibersortx_input==T){
    message('>>> build input for cibersortx')
    if(cibersort_downsample==T){
      downsampled_cells=downsampling(scMeta_train,colnames_of_sample,colnames_of_cellType,scale=cibersort_downsample_scale)
      scMeta_downsampled=scMeta_train[downsampled_cells,]
      sc_raw = scExpr_train[,downsampled_cells]
      colnames(sc_raw)=scMeta_downsampled[,colnames_of_cellType]
      sc_raw=cbind(data.frame(GeneSymbol=rownames(sc_raw)),sc_raw)
      write.table(sc_raw,file = cibersort_input_file_name,sep = '\t',row.names = F)
    }else{
      sc_raw = scExpr_train
      colnames(sc_raw)=scMeta_train[,colnames_of_cellType]
      sc_raw=cbind(data.frame(GeneSymbol=rownames(sc_raw)),sc_raw)
      write.table(sc_raw,file = cibersort_input_file_name,sep = '\t',row.names = F)
    }
  }
  
  message('>>> build heter simulation on testing')
  heter=new.env()
  heter$heter_expr=create_heterSimulation(simulated_frac,scExpr_test,scMeta_test,
                                            colnames_of_cellType,colnames_of_sample,
                                            chunk_size_threshold,
                                            scale_to_million)
  simulation_fold$heter=heter
  
  message('>>> build semiheter simulation on testing')
  semiheter=new.env()
  semiheter$semiheter_expr=create_semiheterSimulation(simulated_frac,scExpr_test,scMeta_test,
                                                      colnames_of_cellType,colnames_of_sample,fixed_cell_type,
                                                      ncells_perSample,chunk_size_threshold_for_fixed_cell)
  simulation_fold$semiheter=semiheter
  
  message('>>> build homo simulation on testing')
  simulation_fold$homo=create_homoSimulation(scExpr_test,scMeta_test,colnames_of_cellType,simulated_frac,ncells_perSample)
  
  return(simulation_fold)
}


########################## deconvolution functions #######################
linseed_result=function(Expr,CellTypeNumber){
  l <- LinseedObject$new(Expr)
  l$calculatePairwiseLinearity()
  l$calculateSpearmanCorrelation()
  l$calculateSignificanceLevel(100)
  # l$significancePlot(0.01)
  l$filterDatasetByPval(0.01)
  # l$svdPlot()
  
  l$setCellTypeNumber(CellTypeNumber)
  l$project("full") # projecting full dataset
  # l$projectionPlot(color="filtered")
  l$project("filtered")
  l$smartSearchCorners(dataset="filtered", error="norm")
  
  l$deconvolveByEndpoints()
  return(l$proportions)
}

CAMfree_result=function(Expr,CellTypeNumber){
  rCAM=CAM(Expr, K = CellTypeNumber)
  rCAM_proportion=Amat(rCAM, CellTypeNumber) %>% t()
  colnames(rCAM_proportion)=colnames(Expr)
  return(rCAM_proportion)
}

gsva_result=function(Expr,gsva_list){
  # filter element in the list with at least 2 elements
  filter_criteria= do.call(c,lapply(gsva_list,length))
  gsva_list=gsva_list[names(filter_criteria)[filter_criteria>=2]]
  gsva_scores=GSVA::gsva(Expr,gsva_list,method='ssgsea',ssgsea.norm=F)
  apply(gsva_scores, 1, scale01) %>% t()# scale the scores to [0,1] for simplicity
}

PCA_result=function(Expr,gsva_list){
  # filter element in the list with at least 2 elements
  filter_criteria= do.call(c,lapply(gsva_list,length))
  gsva_list=gsva_list[names(filter_criteria)[filter_criteria>=2]]
  
  quick_pca=function(genes,Expr){
    mat=Expr[rownames(Expr) %in% genes,]
    exp.PCA<-FactoMineR::PCA(t(as.matrix(mat)),graph = F,scale.unit = T) 
    exp.PCA$ind$coord[,1] 
  }
  m=do.call(rbind,lapply(gsva_list, quick_pca,Expr))
  return(m)
}

avg_result=function(Expr,gsva_list){
  # filter element in the list with at least 2 elements
  filter_criteria= do.call(c,lapply(gsva_list,length))
  gsva_list=gsva_list[names(filter_criteria)[filter_criteria>=2]]
  
  quick_avg=function(genes,Expr){
    mat=Expr[rownames(Expr) %in% genes,]
    colMeans(mat)
  }
  m=do.call(rbind,lapply(gsva_list, quick_avg,Expr))
  return(m)
}

CAMmarker_result=function(Expr,gsva_list){
  est <- tryCatch({CAMTHC::AfromMarkers(Expr,gsva_list,  scaleRecover = TRUE)},
                  error = function(e) NA)
  colnames(est)=names(gsva_list)
  return(t(est))
  # https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/CAMTHC/inst/doc/camthc.html
  # https://github.com/rstudio/rstudio/issues/11076
} 

solveNNLS= function(Y, basis,   lessThanOne = T, weighted=F, rescale=T, quantile=0.9, W=NULL) {
  
  cm=commonRows(Y, basis)
  if(length(cm)<2*ncol(basis)){
    warning("not enough features in intersection")
    stop()
  }
  Y=Y[cm,]
  basis=basis[cm,]
  #rescale to avoid numerical overflow
  if(max(Y,basis)>100){
    ss=100/max(Y,basis)
    Y=ss*Y
    basis=ss*basis
  }  
  if(rescale){
    qqY=quantile(rowMeans(Y), quantile)
    qqB=quantile(rowMeans(basis), quantile)
    basis=basis/qqB*qqY
    
  }
  nCol = dim(basis)[2]
  nSubj = dim(Y)[2]
  mixCoef = matrix(0, nSubj, nCol)
  rownames(mixCoef) = colnames(Y)
  colnames(mixCoef) = colnames(basis)
  if(!weighted){
    Dmat = t(basis) %*% basis
  }
  else if (is.null(W)){
    W=Y
    W[W<0]=0
    W=W+0.01 #need an offset for 0 values
    W=1/W # weights are the inverse of the variance
  }
  
  if (lessThanOne) {
    Amat = cbind(rep(-1, nCol), diag(nCol))
    b0vec = c(-1, rep(0, nCol))
  }
  else {
    Amat = diag(nCol)
    b0vec = rep(0, nCol)
  }
  for (i in 1:nSubj) {
    if (!weighted){
      dvec=t(basis) %*%   Y[,i]
      mixCoef[i, ] = solve.QP(Dmat,dvec, Amat, b0vec)$sol
    }
    else{
      dvec=t(basis) %*%   (Y[,i]*W[,i])
      Dmat = t(basis) %*% (W[,i]*basis)
      mixCoef[i, ] = solve.QP(Dmat,dvec, Amat, b0vec)$sol
    }
  }
  return(mixCoef)
}

deconv_result=function(simulation_fold,Expr,
                       methods=NULL,is_tumor=TRUE,scale_mrna=FALSE, # parameters for immunedeconv::deconvolute
                       run_xcell=F,
                       run_linseed=F,run_CAMfree=F,CellTypeNumber=NULL){

  all_results=list()
 
 if(!is.null(methods)){
    # use methods supported by immunedeconv: https://github.com/omnideconv/immunedeconv
    # example: methods=c('epic','quantiseq','xcell')
    message('>>> Running immunedeconv::deconvolute')
    for (i in 1:length(methods)){
      all_results[[i]]=immunedeconv::deconvolute(Expr,methods[i], column="gene_symbol",scale_mrna = scale_mrna, tumor = is_tumor) %>% column_to_rownames('cell_type') %>% as.matrix()
    }
    names(all_results)=methods
 }
  
  if(run_xcell==T){
    message('>>> Running xCell')
    all_results[['xCellFull']]=xCell::xCellAnalysis(Expr,rnaseq = T) 
  }
  
  message('>>> Running enrichment analysis based on customized gene list')
  gene_list_names=names(simulation_fold$DE_list_bundles)
  for (gene_list_name in gene_list_names){
    all_results[[paste0('MarkerBased_gsva_',gene_list_name)]]=gsva_result(Expr,simulation_fold$DE_list_bundles[[gene_list_name]])
    all_results[[paste0('MarkerBased_eigengene_',gene_list_name)]]=PCA_result(Expr,simulation_fold$DE_list_bundles[[gene_list_name]])
    all_results[[paste0('MarkerBased_CAMTHC_',gene_list_name)]]=CAMmarker_result(Expr,simulation_fold$DE_list_bundles[[gene_list_name]])
    all_results[[paste0('MarkerBased_avgExpr_',gene_list_name)]]=avg_result(Expr,simulation_fold$DE_list_bundles[[gene_list_name]])
  }
  
  message('>>> Running nnls')
  quick_solveNNLS=function(ref_name,Expr,weighted){
    Expr=as.matrix(Expr)
    cms=commonRows(Expr,ref_list[[ref_name]])  %>% as.matrix()
    ref=ref_list[[ref_name]][cms,] %>% as.matrix()
    expr=Expr[cms,]
    message(paste('find',length(cms),'commone genes for',ref_name))
    solveNNLS(expr, ref,weighted = weighted) %>% t()
  }
  
  quick_nnls=function(ref_name,Expr){
    Expr=as.matrix(Expr)
    cms=commonRows(Expr,ref_list[[ref_name]])  %>% as.matrix()
    ref=ref_list[[ref_name]][cms,] %>% as.matrix()
    expr=Expr[cms,]
    nnls_R=do.call(cbind.data.frame,lapply(apply(expr,2,function(x) nnls::nnls(ref,x)), function(y) y$x))
    nnls_R=apply(nnls_R,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(nnls_R) <- colnames(ref)
    nnls_R
  }
  
  ref_names=names(simulation_fold$ref_list)
  ref_list=simulation_fold$ref_list
  for (ref_name in ref_names){
    # all_results[[paste0('RefBased_solveNNLS_',ref_name)]]=quick_solveNNLS(ref_name,Expr,weighted = F)
    all_results[[paste0('RefBased_solveNNLSWeight_',ref_name)]]=quick_solveNNLS(ref_name,Expr,weighted = T) 
    all_results[[paste0('RefBased_nnls_',ref_name)]]=quick_nnls(ref_name,Expr)
  }
  
  if(run_linseed==T){
    message('>>> Running linseed')
    all_results[['linseed']]=linseed_result(Expr,CellTypeNumber)
  }
  
  if(run_CAMfree==T){
    message('>>> Runing CAMfree')
    all_results[['CAMfree']]=CAMfree_result(Expr,CellTypeNumber)
  }
  
  return(all_results)
}

run_deconv<-function(simulation_repeats_obj,
                     methods=NULL,is_tumor=TRUE,scale_mrna=FALSE,
                     run_xcell=F,run_linseed=F,run_CAMfree=F,CellTypeNumber=NULL){
  performance_obj=list()
  for(i in 1:length(simulation_repeats_obj)){
    deconv_fold=list()
    deconv_fold[['homo_deconvResults']]=deconv_result(simulation_repeats_obj[[i]],simulation_repeats_obj[[i]]$homo,methods,is_tumor,scale_mrna,run_xcell,run_linseed,run_CAMfree,CellTypeNumber)
    deconv_fold[['semiheter_deconvResults']]=deconv_result(simulation_repeats_obj[[i]],simulation_repeats_obj[[i]]$semiheter$semiheter_expr,methods,is_tumor,scale_mrna,run_xcell,run_linseed,run_CAMfree,CellTypeNumber)
    deconv_fold[['heter_deconvResults']]=deconv_result(simulation_repeats_obj[[i]],simulation_repeats_obj[[i]]$heter$heter_expr,methods,is_tumor,scale_mrna,run_xcell,run_linseed,run_CAMfree,CellTypeNumber)
    performance_obj[[i]]=deconv_fold
    message(paste('>>>>>>>>>>>>>>>>>>>>>>>>> finish',i,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'))
  }
  return(performance_obj)
}

############### deconvolution performance on simulation repeats object ############
deconv_result<-function(deconv_name,deconv_results,Y,cell_map,mapping_list){
  # E is the estimated fraction, with cell types in columns and samples in rows
  # Y is the true fraction, with cell types in columns and samples in rows
  # cell_map is used for mapping, is built mannually
  # mapping_list is annotated mannually to map cell types in Y and cell types in cell_map
  E=t(deconv_results[[deconv_name]] ) 
  cell_mapSub=cell_map[cell_map$method==deconv_name,]
  m1=gather(Y %>% as.data.frame(),cell_type,true_frac)
  
  stopifnot(all.equal(names(mapping_list),colnames(Y)))
  
  E=as.data.frame(E)
  E$missing_type=NA # add a missing_type column in E, in case that there's no mapping result
  E=as.matrix(E)
  
  # all.equal(names(mapping_list),colnames(Y))  add this checking step later
  
  maxCorName=NULL
  for (i in 1:ncol(Y)){
    # print(paste('evaluate',colnames(Y)[i],'in',deconv_name))
    ee=cell_mapSub[cell_mapSub$ID %in% mapping_list[[i]],]$method_cell_type # names of mapping list must be the same as colnames(y)
    ee=ee[ee%in%colnames(E)]
    if (length(ee)==0){
      maxCorName[i]='missing_type'
    }else if(length(ee)==1){
      maxCorName[i]=ee
    }else{
      maxCorName[i]=colnames(E[,ee])[apply(cor(E[,ee],Y[,i]),2,which.max)]
    }
  }
  
  E=E[,maxCorName]
  colnames(E)=make.unique(colnames(E))
  m2=gather(E %>% as.data.frame(),maxCorName,estimate)
  M=cbind(m1,m2)
  
  summ <- M %>% 
    group_by(cell_type) %>% # note that in summ, cell_type is ordered alphabetically
    summarise(#Rsq = R2(Sepal.Length, Petal.Length),
      RMSE = caret::RMSE(true_frac, estimate,na.rm = T),
      cor = cor(true_frac,estimate,method = 'pearson',use = 'pairwise.complete.obs')) %>% 
    mutate_if(is.numeric, round, digits=2) 
  summ$maxCorName=M$maxCorName[match(summ$cell_type,M$cell_type)]
  return(list(M=M,summ=summ))
}

run_evalu<-function(deconv_repeats_obj,Y,cell_map,mapping_list){
  stopifnot(length(deconv_repeats_obj)>0)
  
  for(i in 1:length(deconv_repeats_obj)){
    d=names(deconv_repeats_obj[[i]][[1]])
    
    get_maxcor_from_R <- function(deconv_name,R){
      R[[deconv_name]]$summ$cor
    }
    
    get_RMSE_from_R <- function(deconv_name,R){
      R[[deconv_name]]$summ$RMSE
    }
    
    listings=names(deconv_repeats_obj[[i]]) 

    for (listing in listings){
      evalu=list()
      evalu[['R']]=lapply(d,deconv_result,deconv_repeats_obj[[i]][[listing]],Y,cell_map,mapping_list)
      names(evalu[['R']])=d
      
      stopifnot(all.equal(evalu[['R']][[1]]$summ$cell_type,colnames(Y)[order(colnames(Y))])) # cell types in summ table should be alphabetically ordered
      
      evalu[['maxcor']]=do.call(cbind,lapply(d,get_maxcor_from_R,evalu[['R']]))
      colnames(evalu[['maxcor']])=d
      rownames(evalu[['maxcor']])=colnames(Y)[order(colnames(Y))]
      evalu[['maxcor']]=as.data.frame(evalu[['maxcor']])
      
      evalu[['RMSE']]=do.call(cbind,lapply(d,get_RMSE_from_R,evalu[['R']]))
      colnames(evalu[['RMSE']])=d
      rownames(evalu[['RMSE']])=colnames(Y)[order(colnames(Y))]
      evalu[['RMSE']]=as.data.frame(evalu[['RMSE']])
      
      evalu_name=paste0(sub("\\_.*", "", listing),'_evalu')
      deconv_repeats_obj[[i]][[evalu_name]]=evalu
      
    }
    
    print(paste('finish',i))
  }
  return(deconv_repeats_obj)
}

quick_gather=function(deconv_performance_repeats,obj_name){
  stopifnot(obj_name %in% c('maxcor','RMSE'))
  all=data.frame(matrix(NA,ncol = 5,nrow=0))
  colnames(all)=c('cell_type','method',obj_name,'group','rep')
  for(i in 1:length(deconv_performance_repeats)){
    homo_long=gather(deconv_performance_repeats[[i]]$homo_evalu[[obj_name]]  %>% rownames_to_column('cell_type'), method, value,-cell_type) %>% mutate(group='homo')
    semiheter_long=gather(deconv_performance_repeats[[i]]$semiheter_evalu[[obj_name]]  %>% rownames_to_column('cell_type'), method, value,-cell_type) %>% mutate(group='semiheter')
    heter_long=gather(deconv_performance_repeats[[i]]$heter_evalu[[obj_name]]  %>% rownames_to_column('cell_type'), method, value,-cell_type) %>% mutate(group='heter')
    combined=rbind(homo_long,semiheter_long) %>% as.data.frame()
    combined=rbind(combined,heter_long) %>% as.data.frame()
    combined$rep=i
    all=rbind(all,combined)
  }
  return(all)
}

################## visualization functions ###########################
deconv_evalu_plot<-function(l,title=NULL,nrow=2){
  M=l$M
  summ=l$summ
  summ=summ[order(summ$cell_type),] # order summ alphabetically
  
  p=ggplot(M,aes(true_frac,estimate))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1,linetype='dotted',color='red')+
    facet_wrap(~cell_type,nrow=nrow)+  # scatter plot will be arranged alphabetically
    theme(aspect.ratio=1)+
    ggtitle(title)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
  
  p + geom_table_npc(data = summ, label = lapply(split(summ, summ$cell_type),
                                                 FUN = function(entry) {subset(entry, select = -cell_type)}),
                     npcx = 0.00, npcy = 1, hjust = 0, vjust = 1, size=3,
                     table.theme = ttheme_gtlight)
}




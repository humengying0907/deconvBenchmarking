
commonRows=function(data1, data2){
  intersect(rownames(data1), rownames(data2))
}

build_ref_matrix<-function(Expr,cell_type_labels){
  stopifnot(ncol(Expr)==length(cell_type_labels))
  group = list()
  for(i in unique(cell_type_labels)){
    group[[i]] <- which(cell_type_labels %in% i)
  }
  C = do.call(cbind, lapply(group,function(x) Matrix::rowMeans(Expr[,x,drop=F])))
  C
}


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

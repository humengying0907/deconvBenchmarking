
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


xnea_fpkm2tpm <- function(fpkmfile,mapfile,to_use_genes,reNormalize = F){
  require(data.table, quietly = T)

  fpkm_expr = fread(fpkmfile) %>% column_to_rownames('Ensembl_ID')
  fpkm_expr = 2^fpkm_expr-1
  fpkmToTpm <- function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  }
  tpm<-apply(fpkm_expr,2,fpkmToTpm)

  gene_mapping = read.delim(mapfile)
  ensembl = gene_mapping$id[match(to_use_genes,gene_mapping$gene)]
  ensembl = ensembl[!is.na(ensembl)]

  tpm = tpm[ensembl,]
  rownames(tpm) = gene_mapping$gene[match(rownames(tpm),gene_mapping$id)]

  if(reNormalize){
    tpm = sweep(tpm,2,colSums(tpm),'/')
  }

  return(tpm)
}

tumor_purity <-function(tcga_barcodes){
  require(TCGAbiolinks, quietly = T)
  data("Tumor.purity")
  tcga_purity=Tumor.purity
  purity_process = function(x){
    x = gsub(',','.',x)
    x = as.numeric(x)
    x[is.na(x)]=NA
    return(x)
  }
  tcga_purity[,3:7]=apply(tcga_purity[,3:7],2,purity_process)

  purity = tcga_purity[tcga_purity$Sample.ID %in% tcga_barcodes,]
  print(paste(nrow(purity),'out of',length(tcga_barcodes),'purity data available'))
  rownames(purity) = NULL
  purity = purity %>% column_to_rownames('Sample.ID')

  purity = append_ABSOLUTE_pancanatlas(purity)
  return(purity)
}

append_ABSOLUTE_pancanatlas = function(tumor_purity_table,ABSOLUTE_pancanatlas = NULL){
  if(is.null(ABSOLUTE_pancanatlas)){
    ABSOLUTE_pancanatlas = read.delim('http://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5')
  }

  tumor_purity_table = as.data.frame(tumor_purity_table)
  find_rowID = function(barcode){
    id = which(grepl(barcode,ABSOLUTE_pancanatlas$sample))
    if(length(id)<1){
      return(NA)
    }else{
      return(id)
    }
  }

  tumor_purity_table$ABSOLUTE_GDC = ABSOLUTE_pancanatlas$purity[do.call(c,lapply(rownames(tumor_purity_table),find_rowID))]
  return(tumor_purity_table)
}


build_tcga_obj = function(tcga_abbreviation,
                          to_use_genes,
                          purity_methods){
  if(is.na(tcga_abbreviation)){
    stop('Please provide the TCGA cohort abbreviation, such as "SKCM", to specifically include this cohort in the bulk object')
  }

  tcga_abbreviation = stringr::str_to_upper(tcga_abbreviation)
  mapfile = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v22.annotation.gene.probeMap'
  fpkmfile = paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-',tcga_abbreviation,'.htseq_fpkm.tsv.gz')

  message('downloading tcga expression from Xena browser')
  tcga_expr = xnea_fpkm2tpm(fpkmfile,mapfile,to_use_genes)
  purity_allMethods = tumor_purity(colnames(tcga_expr))

  purity = matrix(NA,
                  nrow = ncol(tcga_expr),
                  ncol = length(purity_methods),dimnames = list(colnames(tcga_expr),purity_methods))

  for(i in 1:length(purity_methods)){
    purity[,i] = purity_allMethods[,purity_methods[i]][match(colnames(tcga_expr),rownames(purity_allMethods))]
  }

  tcga_obj = list()
  tcga_obj$simulated_bulk = tcga_expr
  tcga_obj$simulated_frac = purity

  return(tcga_obj)
}

get_subcluster = function(scExpr,cell_type_labels, min.subcluster.size = 20){
  require(scran,quietly = T) %>% suppressMessages()

  stopifnot(ncol(scExpr)==length(cell_type_labels))

  group = list()
  for(i in unique(cell_type_labels)){
    group[[i]] <- which(cell_type_labels %in% i)
  }

  get_sublabels = function(ct){
    cl = scran::quickCluster(scExpr[,group[[ct]]],min.size = min(min.subcluster.size,length(group[[ct]])))
    labels = paste0(ct,'_',cl)
    return(labels)
  }

  f = lapply(names(group),get_sublabels)
  names(f) = names(group)

  subcluster_IDtable = data.frame(id = do.call(c,group),
                                  label = do.call(c,f))

  subcluster_IDtable = subcluster_IDtable[order(subcluster_IDtable$id),]
  return(subcluster_IDtable$label)
}


msigdbr_to_list = function(res){
  gs_names = unique(res$gs_name)
  gs_list = list()
  for(gs_name in gs_names){
    gs_list[[gs_name]] = res$gene_symbol[res$gs_name ==gs_name] %>% stringr::str_to_upper()
  }
  return(gs_list)
}

genelist_cv=function(genelist,simulator_Statistics_res){
  summary_stats = simulator_Statistics_res$summary_stats
  summary_stats$gene = stringr::str_to_upper(summary_stats$gene)

  involved_genes=do.call(c,genelist)
  involved_genes=involved_genes[!duplicated(involved_genes)]
  involved_genes=unname(involved_genes)

  df=summary_stats[summary_stats$gene %in% involved_genes,]

  avg_cv_long=do.call(rbind,lapply(genelist,function(x)(df[df$gene %in% x,] %>% group_by(class) %>% summarise(avg_cv=mean(cv,na.rm=T)))))
  avg_cv_long$genelist_name=sub("\\..*","",rownames(avg_cv_long))
  avg_cv_long = data.frame(avg_cv_long)

  return(avg_cv_long)
}




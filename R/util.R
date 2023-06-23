
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

  return(purity)
}

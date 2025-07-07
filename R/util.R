
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


<<<<<<< HEAD

#' Build TCGA bulk expression object for a given cohort
#' @description
#' Downloads TCGA expression data (FPKM) for the specified cohort from UCSC Xena Browser, converts it to TPM, and attaches estimated tumor purity for each sample as an approximation of malignant proportion.
#'
#' @param tcga_abbreviation a character indicating tcga abbreviation for the tcga cohort to include, for example 'SKCM'
#' @param to_use_genes Character vector. Genes to include in the output. If NULL, all available genes are used.
#' @param mapfile_link Character. URL for the gene mapping file (optional; provide only if a custom mapping file is needed).
#' @param fpkmfile_link Character. URL for the FPKM file (optional; provide only if a custom FPKM file is needed).
#'
#' @return a TCGA bulk expression object with TPM expression and purity estimates, which can be used as  an approximation of malignant proportion.
#' @export
#'
#' @examples
#' \dontrun{
#' build_tcga_obj('SKCM')
#' }
build_tcga_obj = function(tcga_abbreviation,
                          to_use_genes = NULL,
                          mapfile_link = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap',
                          fpkmfile_link = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-THCA.star_fpkm.tsv.gz'){
=======
#' Get tcga expression and purity estimates
#' @description Given a tcga_abbreviation, download tcga expression from Xena browser and get purity estimates
#'
#' @param tcga_abbreviation a character indicating tcga abbreviation for the tcga cohort to include, for example 'SKCM'
#' @param to_use_genes a vector indicating which genes to include for the exported tcga expression. Set to NULL to export full expression retrieved from xena browser
#' @param purity_methods a vector indicating tumor purity estimation method that is utilized as a means of estimating the malignant proportion within the exported object for TCGA expression data.
#'     Available methods include 'ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC' and 'CPE' and 'ABSOLUTE_GDC' (ABSOLUTE_GDC contains ABSOLUTE downloaded from https://gdc.cancer.gov/about-data/publications/pancanatlas).
#'    Default = c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE', 'ABSOLUTE_GDC'). Make sure you have the suggested package 'TCGAbiolinks' installed
#' @return a list of tcga expression (named as simulated_bulk, TPM normalized and has hgnc symbols as rownames) and purity estimates (named as simulated_frac)
#' @export
#'
#' @examples
#' \dontrun{build_tcga_obj('SKCM')}
#'
build_tcga_obj = function(tcga_abbreviation,
                          to_use_genes = NULL,
                          purity_methods = c('ESTIMATE', 'ABSOLUTE', 'LUMP', 'IHC', 'CPE','ABSOLUTE_GDC')){
>>>>>>> 7ff7b730a10aa03a7aba7cf0128d54a5c25083f8
  if(is.na(tcga_abbreviation)){
    stop('Please provide the TCGA cohort abbreviation, such as "SKCM", to specifically include this cohort in the bulk object')
  }

  tcga_abbreviation = stringr::str_to_upper(tcga_abbreviation)
  mapfile = 'https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap'
  fpkmfile = paste0('https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-',tcga_abbreviation,'.star_fpkm.tsv.gz')

  map_resp1 <- HEAD(mapfile)

  if(status_code(map_resp1) != 200){
    map_resp2 <- HEAD(mapfile_link)
    if(status_code(map_resp2) != 200){
      stop('Neither the default nor the provided mapfile_link is valida. Please check the UCSC Xena Browser for the most up-to-date link')
    }
  }

  fpkm_resp1 <- httr::HEAD(fpkmfile)
  if (httr::status_code(fpkm_resp1) != 200) {
    fpkm_resp2 <- httr::HEAD(fpkmfile_link)
    if (httr::status_code(fpkm_resp2) != 200) {
      stop('Neither the default nor the provided fpkmfile_link is valida. Please check the UCSC Xena Browser for the most up-to-date link')
    }
  }

  message('downloading tcga expression from Xena browser')
  if(is.null(to_use_genes)){
    gene_mapping = read.delim(mapfile)
    to_use_genes = gene_mapping$gene
  }
  tcga_expr = xnea_fpkm2tpm(fpkmfile,mapfile,to_use_genes)
  purity_allMethods = tumor_purity(colnames(tcga_expr))

  tcga_obj = list()
  tcga_obj$tcga_bulk = tcga_expr
  tcga_obj$tcga_purity = purity_allMethods

  return(tcga_obj)
}

<<<<<<< HEAD

get_subcluster = function(scExpr,cell_type_labels, min.subcluster.size = 20){
  require(scran,quietly = T) %>% suppressMessages()
=======
#' Get sub-clustering labels
#' @description Given single-cell expression and cell-type labels, return sub-clustering labels for each cell type. Currently supported
#'    subclustering-methods include scran::quickCluster() and SCISSORS::ReclusterCells().
#'
#' @param scExpr single cell expression matrix with genes in rows and cells in columns
#' @param cell_type_labels a vector indicating cell-type level annotations
#' @param subcluster_method subclustering methods to use. Available options include 'scran' and 'SCISSORS'. Default = 'scran'.
#' @param min.subcluster.size min.size parameter passed to scran::quickCluster(). Default = 20.
#' @param CalculateSilhouette a logical variable determine whether Silhouette score need to be calculated to determin which cell-types to be reclustered. Default = T. When set to FALSE, will
#'    recluster every cell type present in 'cell_type_labels'
#' @param SilhouetteScores.cutoff a cutoff determining which cell-types to be reclustered with 'SCISSORS' method. Cell-types with SilhouetteScores greater than this
#'    cutoff will not be reclustered.
#' @param ... additional parameters pass to SCISSORS::ReclusterCells()
#'
#' @return a vector of subclustering labels
#' @export
#'
get_subcluster = function(scExpr,cell_type_labels,
                          subcluster_method = 'scran',
                          min.subcluster.size = 20,
                          CalculateSilhouette = T,
                          SilhouetteScores.cutoff = 0.7,...){
>>>>>>> 7ff7b730a10aa03a7aba7cf0128d54a5c25083f8

  cell_type_labels = as.vector(cell_type_labels)

  if (subcluster_method == 'scran'){
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
  }else if(subcluster_method == 'SCISSORS'){
    require(Seurat) %>% suppressMessages()
    require(SCISSORS) %>% suppressMessages()
    seurat_obj = CreateSeuratObject(counts = scExpr)
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = "vst", nfeatures = 5000)
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    seurat_obj <- FindNeighbors(seurat_obj)
    seurat_obj$cell_type = cell_type_labels
    seurat_obj$seurat_clusters = as.factor(as.numeric(as.factor(cell_type_labels)))
    Idents(seurat_obj) = seurat_obj$seurat_clusters
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

    annot = levels(as.factor(cell_type_labels))

    if(length(annot) ==1){
      warning('only one cluster provided, SilhouetteScores is not available, will call ReclusterCells() on the entire scExpr')
      message(paste('recluster',annot))
      reclust_obj  = SCISSORS::ReclusterCells(seurat_obj,
                                              which.clust = 1,
                                              merge.clusters = F,
                                              n.PC = 15,
                                              use.parallel = FALSE,
                                              redo.embedding = T,...)
      subcluster_label = paste0(annot,'_',reclust_obj$seurat_clusters)
      gc()
      return(subcluster_label)
    }else{
      if(CalculateSilhouette){
        # decide which clusters to recluster
        SilhouetteScores = ComputeSilhouetteScores(seurat_obj,avg = T) %>% as.numeric()
        recluster_id = which(SilhouetteScores < SilhouetteScores.cutoff)

        gc()

        if(length(recluster_id)<1){
          message('no cell-type need to be reclustered under the given SilhouetteScores.cutoff, will return original cell type labels')
          return(cell_type_labels)
        }

      }else{
        recluster_id = seq(1,length(unique(cell_type_labels)))
      }

      seurat_obj$subcluster = seurat_obj$cell_type

      for(i in recluster_id){
        message(paste('recluster',annot[i]))
        reclust_obj  = SCISSORS::ReclusterCells(seurat_obj,
                                                which.clust = i,
                                                merge.clusters = F,
                                                n.PC = 15,
                                                use.parallel = FALSE,
                                                redo.embedding = T,...)
        id = match(colnames(reclust_obj),colnames(seurat_obj))
        subcluster_label = paste0(annot[i],'_',reclust_obj$seurat_clusters)

        seurat_obj$subcluster[id] = subcluster_label
        gc()

      }
      return(seurat_obj$subcluster)
    }
  }else{
    stop('please provide valid subcluster method')
  }
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



